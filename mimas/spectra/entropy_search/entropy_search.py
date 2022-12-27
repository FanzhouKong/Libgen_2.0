import numpy as np
import struct
from mimas.spectra.spectral_entropy_c.functions import apply_weight_to_intensity
from mimas.helper.multiplecore import convert_numpy_array_to_shared_memory

from .entropy_search_core_fast import entropy_similarity_search_fast


class EntropySearchCore:
    """
    This class will only do the core function of entropy search.
    This class only calculate the product ions
    """

    def __init__(self):
        self.total_spectra_num = 0
        self.index = []

    def build_index(self, all_peaks_list):
        """
        Build the index for the spectra library.
        :param all_peaks_list: list of all peaks in the spectra library. 
                                Each element is a 2-D numpy array of shape (n, 2).
                                The n is the number of peaks in the spectra, n >= 1.
                                The sum of all intensities in the spectra is 1.
        """
        if len(all_peaks_list) == 0:
            self.total_spectra_num = 0
            self.index = []
            return self.index
        all_entropy, all_half_entropy = [], []
        all_peaks_info = []
        for idx, peaks in enumerate(all_peaks_list):
            assert peaks.shape[0] >= 1  # Check if the number of peaks is greater than 0.
            assert np.sum(peaks[:, 1]) >= 0.999  # Check if the sum of all intensities is 1.
            assert np.all(peaks[:-1, 0] <= peaks[1:, 0])  # Check if the peaks are sorted.

            # Pre-calculate library peaks entropy
            peaks_clean = np.asarray(apply_weight_to_intensity(peaks))
            peaks_spectral_entropy = -np.sum(peaks_clean[:, 1]*np.log2(peaks_clean[:, 1]))
            peaks_clean[:, 1] /= 2
            peaks_half_spectral_entropy = -np.sum(peaks_clean[:, 1]*np.log2(peaks_clean[:, 1]))

            # Store the precursor m/z and the spectral entropy
            all_entropy.append(peaks_spectral_entropy)
            all_half_entropy.append(peaks_half_spectral_entropy)

            # Store the peaks
            peaks_new = np.concatenate((peaks_clean, np.zeros((peaks_clean.shape[0], 1), dtype=np.int64)+idx), axis=1, dtype=object)
            all_peaks_info.append(peaks_new)

        # Generate the index
        all_peaks_info = np.concatenate(all_peaks_info, axis=0)
        all_peaks_info = all_peaks_info[all_peaks_info[:, 0].argsort()]

        all_peaks_mz = np.floor(all_peaks_info[:, 0]*1000).astype(np.int64)
        all_peaks_mz_idx_all = np.array(range(int(np.max(all_peaks_mz))), dtype=np.int64)
        all_peaks_mz_idx_start = np.concatenate([np.searchsorted(all_peaks_mz, all_peaks_mz_idx_all), [len(all_peaks_info)]], dtype=np.int64).astype(np.int64)

        all_peaks_mz = all_peaks_info[:, 0].astype(np.float32)
        all_peaks_intensity = all_peaks_info[:, 1].astype(np.float32)
        all_peaks_spec_idx = all_peaks_info[:, 2].astype(np.int64)

        all_entropy, all_half_entropy = np.array(all_entropy, dtype=np.float32), np.array(all_half_entropy, dtype=np.float32)

        self.total_spectra_num = len(all_entropy)
        self.index = all_entropy, all_half_entropy, all_peaks_mz_idx_start, all_peaks_mz, all_peaks_intensity, all_peaks_spec_idx
        return self.index

    def search(self, *, peaks, ms2_tolerance_in_da,
               search_type=0, search_spectra_idx_min=0, search_spectra_idx_max=0, search_array=None):
        """
        search_type is 0: search all spectra.
        search_type is 1: search spectra in the range [search_spectra_idx_min, search_spectra_idx_max).
        search_type is 2: search spectra in the array search_array with entry equals 1, the length of search_array should be equal to the self.total_spectra_num
        """
        if self.total_spectra_num == 0:
            return np.zeros(0, dtype=np.float32)
        # Calculate query peaks entropy
        peaks = np.asarray(apply_weight_to_intensity(peaks))
        query_entropy = -np.sum(peaks[:, 1]*np.log2(peaks[:, 1]))
        peaks[:, 1] /= 2
        query_half_entropy = -np.sum(peaks[:, 1]*np.log2(peaks[:, 1]))
        mixed_spectra_entropy = np.zeros(self.total_spectra_num, dtype=np.float32)

        library_entropy, library_half_entropy, library_peaks_mz_idx_start, library_peaks_mz, library_peaks_intensity, library_peaks_spec_idx = self.index

        # Go through all the peaks in the spectrum
        for mz, intensity_half in peaks:
            # Determine the mz index range
            product_mz_idx_min = calculate_index_for_product_ion_mz(mz - ms2_tolerance_in_da, library_peaks_mz, library_peaks_mz_idx_start, 'left')
            product_mz_idx_max = calculate_index_for_product_ion_mz(mz + ms2_tolerance_in_da, library_peaks_mz, library_peaks_mz_idx_start, 'right')

            entropy_similarity_search_fast(product_mz_idx_min, product_mz_idx_max, intensity_half, mixed_spectra_entropy,
                                           library_peaks_intensity, library_peaks_spec_idx,
                                           search_type, search_spectra_idx_min, search_spectra_idx_max, search_array)

        entropy_spectra_ab = library_half_entropy + query_half_entropy + mixed_spectra_entropy
        entropy_similarity = 1 - (entropy_spectra_ab - (library_entropy + query_entropy) / 2)
        return entropy_similarity

    def get_total_spectra_num(self):
        return self.total_spectra_num

    def write_to_file(self, fo):
        library_entropy, library_half_entropy, library_peaks_mz_idx_start, library_peaks_mz, library_peaks_intensity, library_peaks_spec_idx = self.index
        len_library_peaks_mz_idx_start = len(library_peaks_mz_idx_start)
        total_peaks_num = len(library_peaks_mz)
        fo.write(struct.pack('QQQ', self.total_spectra_num, total_peaks_num,  len_library_peaks_mz_idx_start))
        for data in self.index:
            data.tofile(fo)

    def read_from_file(self, fi):
        self.total_spectra_num, total_peaks_num, len_library_peaks_mz_idx_start = struct.unpack('QQQ', fi.read(8+8))
        library_entropy = np.fromfile(fi, dtype=np.float32, count=self.total_spectra_num)
        library_half_entropy = np.fromfile(fi, dtype=np.float32, count=self.total_spectra_num)
        library_peaks_mz_idx_start = np.fromfile(fi, dtype=np.int64, count=len_library_peaks_mz_idx_start)
        library_peaks_mz = np.fromfile(fi, dtype=np.float32, count=total_peaks_num)
        library_peaks_intensity = np.fromfile(fi, dtype=np.float32, count=total_peaks_num)
        library_peaks_spec_idx = np.fromfile(fi, dtype=np.int64, count=total_peaks_num)
        self.index = library_entropy, library_half_entropy, library_peaks_mz_idx_start, library_peaks_mz, library_peaks_intensity, library_peaks_spec_idx

    def move_index_array_to_shared_memory(self):
        if len(self.index) == 0:
            return

        library_entropy, library_half_entropy, library_peaks_mz_idx_start, library_peaks_mz, library_peaks_intensity, library_peaks_spec_idx = self.index

        library_entropy, library_half_entropy, library_peaks_mz, library_peaks_intensity = \
            convert_numpy_array_to_shared_memory(library_entropy, "f"),\
            convert_numpy_array_to_shared_memory(library_half_entropy, "f"),\
            convert_numpy_array_to_shared_memory(library_peaks_mz, "f"),\
            convert_numpy_array_to_shared_memory(library_peaks_intensity, "f")
        library_peaks_mz_idx_start, library_peaks_spec_idx = \
            convert_numpy_array_to_shared_memory(library_peaks_mz_idx_start, "q"),\
            convert_numpy_array_to_shared_memory(library_peaks_spec_idx, "q")

        self.index = library_entropy, library_half_entropy, library_peaks_mz_idx_start, library_peaks_mz, library_peaks_intensity, library_peaks_spec_idx


def calculate_index_for_product_ion_mz(product_mz, library_peaks_mz, library_peaks_mz_idx_start, side):
    product_mz_min_int = (np.floor(product_mz*1000)).astype(int)
    product_mz_max_int = product_mz_min_int + 1

    if product_mz_min_int >= len(library_peaks_mz_idx_start):
        product_mz_idx_search_start = library_peaks_mz_idx_start[-1].astype(int)
    else:
        product_mz_idx_search_start = library_peaks_mz_idx_start[product_mz_min_int].astype(int)

    if product_mz_max_int >= len(library_peaks_mz_idx_start):
        product_mz_idx_search_end = library_peaks_mz_idx_start[-1].astype(int)
    else:
        product_mz_idx_search_end = library_peaks_mz_idx_start[product_mz_max_int].astype(int)

    return product_mz_idx_search_start + np.searchsorted(
        library_peaks_mz[product_mz_idx_search_start:product_mz_idx_search_end], product_mz, side=side)


def entropy_similarity_search(product_mz_idx_min, product_mz_idx_max, intensity, mixed_spectra_entropy,
                              library_peaks_intensity, library_spec_idx_array,
                              search_type, search_spectra_idx_min, search_spectra_idx_max, search_array):
    """
    The mixed_spectra_entropy will be modified in this function.
    search_type is 0: search all spectra.
    search_type is 1: search spectra in the range [search_spectra_idx_min, search_spectra_idx_max).
    search_type is 2: search spectra in the array search_array with entry equals 1, the length of search_array should be equal to the self.total_spectra_num

    Note: the intensity here should be half of the original intensity.
    """
    for idx in range(product_mz_idx_min, product_mz_idx_max):
        library_spec_idx = library_spec_idx_array[idx]
        if (search_type == 0) or \
            (search_type == 1 and search_spectra_idx_min <= library_spec_idx and library_spec_idx < search_spectra_idx_max) or \
                (search_type == 2 and search_array[library_spec_idx]):
            # Match this peak
            library_peak_intensity = library_peaks_intensity[idx]
            intensity_ab = intensity + library_peak_intensity

            mixed_spectra_entropy[library_spec_idx] -= \
                intensity_ab * np.log2(intensity_ab) - \
                intensity * np.log2(intensity) - \
                library_peak_intensity * np.log2(library_peak_intensity)


class EntropySearch(EntropySearchCore):
    """
    This class consider the precursor ion and product ions
    """

    def __init__(self):
        super().__init__()
        self.precursor_mz: np.ndarray = np.zeros(0, dtype=np.float32)

    def build_index(self, all_spectra, sort_by_precursor_mz=True):
        """
        This function will build the index from the spectra generator.
        The spectra need to be sorted by precursor m/z.
        The peaks need to be cleaned and normalize to 1 before using this function.
        This function will use the input as is, and will not do any pre-processing except the intensity weighting by spectral entropy.
        The following keys are used:
            "precursor_mz": precursor m/z
            "peaks": a numpy array of the peaks, with the first column as the m/z and the second column as the intensity

        """
        if sort_by_precursor_mz:
            all_spectra_ori = all_spectra
            all_spectra = [x for x in all_spectra_ori]
            all_spectra.sort(key=lambda x: x["precursor_mz"])

        self.precursor_mz = np.array([x["precursor_mz"] for x in all_spectra], dtype=np.float32)
        super().build_index([x["peaks"] for x in all_spectra])
        return all_spectra

    def search_identity(self, *, precursor_mz, peaks, ms1_tolerance_in_da, ms2_tolerance_in_da):
        # Determine the precursor m/z range
        precursor_mz_min = precursor_mz - ms1_tolerance_in_da
        precursor_mz_max = precursor_mz + ms1_tolerance_in_da
        spectra_idx_min = np.searchsorted(self.precursor_mz, precursor_mz_min, side='left')
        spectra_idx_max = np.searchsorted(self.precursor_mz, precursor_mz_max, side='right')
        if spectra_idx_min >= spectra_idx_max:
            return np.zeros(self.total_spectra_num, dtype=np.float32)
        else:
            entropy_similarity = self.search(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da,
                                             search_type=1, search_spectra_idx_min=spectra_idx_min, search_spectra_idx_max=spectra_idx_max)
            return entropy_similarity

    def search_open(self, *, peaks, ms2_tolerance_in_da):
        entropy_similarity = self.search(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da,
                                         search_type=0)
        return entropy_similarity

    def move_index_array_to_shared_memory(self):
        super().move_index_array_to_shared_memory()

        if len(self.precursor_mz) == 0:
            return

        self.precursor_mz = convert_numpy_array_to_shared_memory(self.precursor_mz, "f")
