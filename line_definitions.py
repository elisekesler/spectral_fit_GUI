from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

@dataclass
class LineDefinition:
    name: str  # Display name (e.g., "Lya")
    rest_wavelength: float  # Rest wavelength in Angstroms
    formal_name: str  # More formal name (e.g., "Lyman-alpha")
    is_doublet: bool = False
    secondary_wavelength: Optional[float] = None
    doublet_ratio: Optional[float] = 0.5
    category: str = "other"  # UV, optical, etc.

# UV lines
UV_LINES = [
    LineDefinition("Lya", 1215.67, "Lyman-alpha", category="UV"),
    LineDefinition("OVI", 1031.92, "Oxygen VI", True, 1037.61, 0.5, category="UV"),
    LineDefinition("NV", 1238.82, "Nitrogen V", True, 1242.8, 0.5, category="UV"),
    LineDefinition("SiIV", 1393.75, "Silicon IV", True, 1402.77, 0.5, category="UV"),
    LineDefinition("CIV", 1548.19, "Carbon IV", True, 1550.77, 0.5, category="UV"),
    LineDefinition("HeII", 1640.42, "Helium II", category="UV"),
]

# Optical lines
OPTICAL_LINES = [
    LineDefinition("Ha", 6562.8, "Hydrogen-alpha", category="optical"),
    LineDefinition("Hb", 4861.33, "Hydrogen-beta", category="optical"),
    LineDefinition("[OIII]", 4958.91, "Oxygen III", True, 5006.84, (1/2.95), category="optical"),
    LineDefinition("[NII]", 6548.05, "Nitrogen II", True, 6583.45, (1/2.95), category="optical"),
    LineDefinition("HeI", 5875.67, "Helium I", category="optical"),
]

# All lines
ALL_LINES = UV_LINES + OPTICAL_LINES

# Map names to definitions for quick lookup
LINE_MAP = {line.name: line for line in ALL_LINES}