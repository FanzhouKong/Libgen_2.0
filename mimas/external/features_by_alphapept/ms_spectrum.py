import numpy as np
import struct


class MSSpectrum:
    def __init__(self) -> None:
        self.scan: int
        self.precursor_mz: float
        self.charge: int
        self.peaks: np.ndarray

    def from_input(self, scan: int, precursor_mz: float, charge: int, peaks: np.ndarray) -> None:
        self.scan = scan
        self.precursor_mz = precursor_mz
        self.charge = charge
        self.peaks = peaks

    def read_from_file(self, f):
        self.scan, self.precursor_mz,  self.charge, peak_len = struct.unpack(r'IfhH', f.read(12))
        self.peaks = np.fromfile(f, dtype=np.float32, count=2*peak_len).reshape((-1, 2))
        pass

    def write_to_file(self, f):
        info = struct.pack(r'IfhH', self.scan, self.precursor_mz,  self.charge, self.peaks.shape[0])
        f.write(info)
        self.peaks.astype(np.float32).tofile(f)
        return
