class Distances:
    """Container for all distance information that has been extracted from a
    system.
    """

    def __init__(
        self,
        disp_tensor,
        disp_factors,
        dist_matrix,
        dist_matrix_radii,
        cell_list=None,
    ):
        self.disp_tensor = disp_tensor
        self.disp_factors = disp_factors
        self.dist_matrix = dist_matrix
        self.dist_matrix_radii = dist_matrix_radii
        self.cell_list = cell_list
