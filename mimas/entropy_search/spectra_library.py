import numpy as np
import pickle


class SpectraLibrary:
    def __init__(self, parameter):
        """Exmaple of the parameter:
        {
            "calibrate_mz": True,
            "clean_spectra": True,
            "ms2_tolerance_in_da": 0.02,
            "simplify_adduct": False,
            "index_for_hybrid_search": True,
        }
        """
        # Try to load index file
        self.name = ""
        self.library = {}
        self.metadata = []
        self.parameter = parameter

    def build_index_from_spectral_library(self,  name, file_spectra):
        # Those packages are used only when build index.
        from mimas.file_io import spec_file

        self.name = name
        self.spectra = {1: [], -1: []}
        collections_information = []

        for spec_info in spec_file.read_one_spectrum(file_spectra, ms2_only=True):
            try:
                spec_file.standardize_spectrum(spec_info, {
                    "precursor_mz": [["precursormz"], -1, float],
                    "adduct": [["precursortype", "precursor_type"], "", str],
                    "name": [[],  "", str],
                    "id": [["spectrum_id", "db#", "spectrumid"], "", str],
                    "inchikey": [[], None, None],
                    "smiles": [[], "", str],
                })
            except Exception as e:
                print("Error in parse: ", spec_info)
                continue

            # Organize spectra data
            self._index_one_spec(spec_info, self.spectra, collections_information)

        # Organize collections
        for charge in self.spectra:
            collection = self.spectra[charge]
            mz = np.array([x[0] for x in collection], dtype=np.float32)
            mz_idx = np.argsort(mz)
            collection = [collection[i] for i in mz_idx]
            self.library[charge] = [mz[mz_idx], collection]
        if self.parameter["index_for_hybrid_search"]:
            self.add_index_for_hybrid_search(min_intensity=0.05, max_peak_number=10)
        pass

    def _index_one_spec(self, spec_info: dict, collections_spec: dict, collections_information):
        from mimas.spectra.similarity import clean_spectrum
        from mimas.chem.adduct import Adduct
        from rdkit.Chem import Descriptors
        from rdkit import Chem
        precursor_mz = spec_info["precursor_mz"]
        adduct = spec_info["adduct"]
        peaks = [x for x in spec_info["peaks"] if x[0] > 0 and x[1] > 0]
        peaks = clean_spectrum(peaks,
                               max_mz=precursor_mz-1.6, noise_threshold=0.01, remove_isotope=True,
                               normalize_intensity=True,
                               ms2_da=self.parameter["ms2_tolerance_in_da"])

        # Only process MS2 spectra with: precursor_mz, adduct, peaks.
        if (len(peaks) > 0) and \
                (len(adduct) > 1) and (adduct[-1] in {"+", "-"}) and \
                (precursor_mz > 0):
            adduct_info = Adduct(adduct)
            charge = adduct_info.charge
            if charge not in self.spectra:
                self.spectra[charge] = []

            if self.parameter["calibrate_mz"] and spec_info["smiles"] and adduct_info.correct:
                try:
                    precursor_mz_calibrated = -1
                    mol = Chem.MolFromSmiles(spec_info["smiles"])
                    if mol is not None:
                        mass = Descriptors.ExactMolWt(mol)
                        precursor_mz_calibrated = adduct_info.get_mz_with_adduct(mass)
                        if abs(precursor_mz_calibrated - precursor_mz) < 0.5:
                            precursor_mz = precursor_mz_calibrated
                        else:
                            print(f"Precursor m/z is not calibrated: {precursor_mz} -> {precursor_mz_calibrated} in {spec_info['id']} with adduct {adduct}")
                except:
                    print(f"Precursor m/z is not calibrated: {precursor_mz} -> {precursor_mz_calibrated} in {spec_info['id']} with adduct {adduct}")
                    pass
            if not adduct_info.correct:
                print(f"Adduct is not correct: {adduct} in {spec_info['id']}")

            self.spectra[charge].append((precursor_mz, peaks, len(self.metadata)))

            # Remove adduct mass from the precursor_mz
            if self.parameter["simplify_adduct"] and adduct_info.correct:
                if abs(charge) == 1 and (adduct not in {"[M+H]+", "[M-H]-"}):
                    precursor_mz_simple_adduct = adduct_info.simplify_adduct(mz=precursor_mz)
                    if abs(precursor_mz_simple_adduct - precursor_mz) > 0.5:
                        self.spectra[charge].append((precursor_mz_simple_adduct, peaks, len(self.metadata)))

            info_metadata = {
                "id": spec_info["id"],
                "precursor_mz": precursor_mz,
                "adduct": adduct,
                "name": spec_info["name"],
                "inchikey": spec_info["inchikey"],
                "smiles": spec_info["smiles"],
            }
            for key in {"_scan_number", "_ms_level", "inchikey", "id", "adduct", "name", "smiles", "num peaks", "precursor_mz", "peaks", "inchi"}:
                if key in spec_info:
                    spec_info.pop(key)
            info_metadata["metadata"] = spec_info
            self.metadata.append(info_metadata)
        pass

    def read_index_from_file(self, file_index):
        data = pickle.load(open(file_index, "rb"))
        self.name, self.parameter, self.library, self.metadata = data

    def write_index_to_file(self, file_index):
        data = [self.name, self.parameter, self.library, self.metadata]
        pickle.dump(data, open(file_index, "wb"))

    def get_candidate_spectra_all(self, charge: int):
        if charge in self.library:
            all_spec = self.library[charge][1]
        else:
            all_spec = []
        return all_spec

    def get_candidate_spectra(
            self,
            precursor_mz: float,
            charge: int,
            delta_ppm: float = None,
            delta_da: float = None):
        if delta_ppm:
            delta_da = precursor_mz * delta_ppm * 1e-6

        mz_min = precursor_mz - delta_da
        mz_max = precursor_mz + delta_da

        precursor_mz_idx, spectra = self.library[charge][0:2]
        idx_left = np.searchsorted(precursor_mz_idx, mz_min, side="left")
        idx_right = np.searchsorted(precursor_mz_idx, mz_max, side="right")

        return spectra[idx_left: idx_right]

    def get_metadata(self, library_id):
        return self.metadata[library_id]

    def get_library_name(self):
        return self.name

    def add_index_for_hybrid_search(self, min_intensity=0.05, max_peak_number=10):
        ms2_da = self.parameter["ms2_tolerance_in_da"]

        def convert_idx_to_numpy_array(idx_list):
            idx_list = [np.array(list(x), dtype=np.uint32) for x in idx_list]
            idx_list_len = np.array([len(x) for x in idx_list])
            idx_list_len_cumsum = np.concatenate([np.array([0]), np.cumsum(idx_list_len)], dtype=np.int64)
            idx_list_array = np.concatenate(idx_list)
            return idx_list_len_cumsum, idx_list_array
        for charge in self.library:
            if abs(charge) == 1:
                all_spec = self.library[charge][1]
                all_mz_idx = []
                all_neutral_loss_idx = []
                for i, (precursor_mz, spec, _) in enumerate(all_spec):
                    if len(spec) == 0:
                        continue

                    spec_idx = spec[spec[:, 1] > min_intensity]
                    if len(spec_idx) > max_peak_number:
                        spec_idx = spec_idx[np.argsort(-spec_idx[:, 1])][:max_peak_number]

                    # Add index for m/z
                    mz_idx = spec_idx[:, 0] // ms2_da
                    for idx in mz_idx:
                        idx = int(idx)
                        if idx >= len(all_mz_idx):
                            all_mz_idx += [set() for _ in range(idx - len(all_mz_idx) + 1)]
                        all_mz_idx[idx].add(i)

                    # Add index for neutral loss
                    nl_idx = (precursor_mz - spec_idx[:, 0]) // ms2_da
                    for idx in nl_idx:
                        idx = int(idx)
                        if idx >= len(all_neutral_loss_idx):
                            all_neutral_loss_idx += [set() for _ in range(idx - len(all_neutral_loss_idx) + 1)]
                        all_neutral_loss_idx[idx].add(i)

                all_mz_idx_info, all_mz_idx_data = convert_idx_to_numpy_array(all_mz_idx)
                all_neutral_loss_info, all_neutral_loss_idx = convert_idx_to_numpy_array(all_neutral_loss_idx)
                self.library[charge] = [self.library[charge][0], self.library[charge][1],
                                        ((all_mz_idx_info, all_mz_idx_data), (all_neutral_loss_info, all_neutral_loss_idx))]
            else:
                self.library[charge] = [self.library[charge][0], self.library[charge][1], ([], [])]

    def get_candidate_spectra_for_hybrid_search(self,
                                                precursor_mz: float,
                                                charge: int,
                                                peaks: np.ndarray,
                                                min_intensity=0.05, max_peak_number=10) -> list:
        if len(peaks) == 0:
            return []
        ms2_da = self.parameter["ms2_tolerance_in_da"]
        precursor_mz_idx, all_spec, (library_mz_idx, library_nl_idx) = self.library[charge][0:3]
        result = list()

        peaks_idx = peaks[peaks[:, 1] > min_intensity]
        if len(peaks_idx) > max_peak_number:
            peaks_idx = peaks_idx[np.argsort(-peaks_idx[:, 1])][:max_peak_number]
        for mz in peaks_idx[:, 0]:
            mz_idx = int(mz // ms2_da)
            try:
                range = library_mz_idx[0][(mz_idx-1):(mz_idx+3)]
                if len(range) > 1:
                    result.append(library_mz_idx[1][range[0]:range[-1]])
            except Exception as e:
                pass

        peaks_idx = precursor_mz - peaks_idx[:, 0]
        for nl in peaks_idx:
            nl_idx = int(nl // ms2_da)
            try:
                range = library_nl_idx[0][(mz_idx-1):(mz_idx+3)]
                if len(range) > 1:
                    result.append(library_nl_idx[1][range[0]:range[-1]])
            except Exception as e:
                pass

        result = np.unique(np.concatenate(result))
        return [all_spec[x] for x in result]

    def pop_metadata(self):
        metadata = self.metadata
        self.set_metadata([{} for _ in range(len(metadata))])
        return metadata

    def set_metadata(self, metadata):
        self.metadata = metadata


if __name__ == '__main__':
    from pathlib import Path
    # Test
    filename_spec_1 = r"result/2022_05/0506_hybrid_search_on_nist20_linear_spectra/data/spec_test_1.msp"
    filename_spec_1_index_spectra = Path(filename_spec_1).with_suffix(".spectra.index")
    filename_spec_1_index_metadata = Path(filename_spec_1).with_suffix(".metadata.index")

    db = SpectraLibrary(parameter={
        "calibrate_mz": True,
        "clean_spectra": True,
        "ms2_tolerance_in_da": 0.02,
        "index_for_hybrid_search": True,
    })
    db.build_index_from_spectral_library(name="test", file_spectra=filename_spec_1)
    db.write_spectra_to_file(file_spectra=filename_spec_1_index_spectra)
    db.write_metadata_to_file(file_metadata=filename_spec_1_index_metadata)

    db = SpectraLibrary(parameter={
        "calibrate_mz": True,
        "clean_spectra": True,
        "ms2_tolerance_in_da": 0.02,
        "index_for_hybrid_search": True,
    })
    db.read_spectra_from_file(file_spectra=filename_spec_1_index_spectra)
    db.read_metadata_from_file(file_metadata=filename_spec_1_index_metadata)
    spec = db.get_candidate_spectra(precursor_mz=263.1188, charge=-1, delta_da=3)
    print(db.get_metadata(library_id=spec[0][2]))
