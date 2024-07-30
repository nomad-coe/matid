"use strict";(self.webpackChunkdocs=self.webpackChunkdocs||[]).push([[967],{21:(e,n,s)=>{s.r(n),s.d(n,{assets:()=>c,contentTitle:()=>l,default:()=>d,frontMatter:()=>r,metadata:()=>a,toc:()=>o});var t=s(5893),i=s(1151);const r={sidebar_position:3,sidebar_label:"SymmetryAnalyzer Class"},l="SymmetryAnalyzer Class",a={id:"reference/symmetryanalyzer",title:"SymmetryAnalyzer Class",description:"A base class for getting symmetry related properties of unit cells.",source:"@site/docs/reference/symmetryanalyzer.md",sourceDirName:"reference",slug:"/reference/symmetryanalyzer",permalink:"/matid/docs/reference/symmetryanalyzer",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:3,frontMatter:{sidebar_position:3,sidebar_label:"SymmetryAnalyzer Class"},sidebar:"docsSidebar",previous:{title:"Cluster Class",permalink:"/matid/docs/reference/cluster"},next:{title:"Geometry Module",permalink:"/matid/docs/reference/geometry"}},c={},o=[{value:"__init__",id:"__init__",level:2},{value:"set_system",id:"set_system",level:2},{value:"reset",id:"reset",level:2},{value:"get_material_id",id:"get_material_id",level:2},{value:"get_space_group_number",id:"get_space_group_number",level:2},{value:"get_space_group_international_short",id:"get_space_group_international_short",level:2},{value:"get_hall_symbol",id:"get_hall_symbol",level:2},{value:"get_hall_number",id:"get_hall_number",level:2},{value:"get_point_group",id:"get_point_group",level:2},{value:"get_is_chiral",id:"get_is_chiral",level:2},{value:"get_has_free_wyckoff_parameters",id:"get_has_free_wyckoff_parameters",level:2},{value:"get_crystal_system",id:"get_crystal_system",level:2},{value:"get_bravais_lattice",id:"get_bravais_lattice",level:2},{value:"get_primitive_system",id:"get_primitive_system",level:2},{value:"get_conventional_system",id:"get_conventional_system",level:2},{value:"get_rotations",id:"get_rotations",level:2},{value:"get_translations",id:"get_translations",level:2},{value:"get_choice",id:"get_choice",level:2},{value:"get_wyckoff_letters_original",id:"get_wyckoff_letters_original",level:2},{value:"get_equivalent_atoms_original",id:"get_equivalent_atoms_original",level:2},{value:"get_wyckoff_letters_conventional",id:"get_wyckoff_letters_conventional",level:2},{value:"get_wyckoff_sets_conventional",id:"get_wyckoff_sets_conventional",level:2},{value:"get_equivalent_atoms_conventional",id:"get_equivalent_atoms_conventional",level:2},{value:"get_wyckoff_letters_primitive",id:"get_wyckoff_letters_primitive",level:2},{value:"get_equivalent_atoms_primitive",id:"get_equivalent_atoms_primitive",level:2},{value:"get_symmetry_dataset",id:"get_symmetry_dataset",level:2},{value:"get_symmetry_operations",id:"get_symmetry_operations",level:2}];function h(e){const n={code:"code",em:"em",h1:"h1",h2:"h2",li:"li",p:"p",pre:"pre",strong:"strong",ul:"ul",...(0,i.a)(),...e.components};return(0,t.jsxs)(t.Fragment,{children:[(0,t.jsx)(n.h1,{id:"symmetryanalyzer-class",children:"SymmetryAnalyzer Class"}),"\n",(0,t.jsx)(n.p,{children:"A base class for getting symmetry related properties of unit cells."}),"\n",(0,t.jsx)(n.h2,{id:"__init__",children:"__init__"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def __init__(system=None, symmetry_tol=None, min_2d_thickness=1)\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Arguments"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"system(ASE.Atoms)"})," - The system to inspect."]}),"\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"symmetry_tol(float)"})," - The tolerance for the symmetry detection."]}),"\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"min_2d_thickness(float)"})," - The minimum thickness in angstroms for the\nconventional cell that is returned for 2D systems."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"set_system",children:"set_system"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def set_system(system)\n"})}),"\n",(0,t.jsx)(n.p,{children:"Sets a new system for analysis."}),"\n",(0,t.jsx)(n.h2,{id:"reset",children:"reset"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def reset()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Used to reset all the cached values."}),"\n",(0,t.jsx)(n.h2,{id:"get_material_id",children:"get_material_id"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_material_id()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Returns a 28-character identifier for this material. The identifier\nis calculated by hashing a set of the symmetry properties found in the\nmaterial, including:"}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsx)(n.li,{children:"Space group number"}),"\n",(0,t.jsx)(n.li,{children:"Wyckoff position letters and the species occupied in them"}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_space_group_number",children:"get_space_group_number"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_space_group_number()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"int"})," - The space group number."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_space_group_international_short",children:"get_space_group_international_short"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_space_group_international_short()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"str"})," - The international space group short symbol."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_hall_symbol",children:"get_hall_symbol"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_hall_symbol()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"str"})," - The Hall symbol."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_hall_number",children:"get_hall_number"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_hall_number()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"int"})," - The Hall number."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_point_group",children:"get_point_group"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_point_group()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Symbol of the crystallographic point group in the Hermann-Mauguin\nnotation."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"str"})," - point group symbol"]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_is_chiral",children:"get_is_chiral"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_is_chiral()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Returns a boolean value that tells if this object is chiral or not\n(achiral). A chiral object has symmetry operations that are all proper,\ni.e. their determinant is +1."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"bool"})," - is the object chiral."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_has_free_wyckoff_parameters",children:"get_has_free_wyckoff_parameters"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_has_free_wyckoff_parameters()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Tells whether this system has Wyckoff positions with free variables."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"bool"})," - Indicates the presence of Wyckoff positions with free variables."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_crystal_system",children:"get_crystal_system"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_crystal_system()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get the crystal system based on the space group number. There are\nseven different crystal systems:"}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsx)(n.li,{children:"Triclinic"}),"\n",(0,t.jsx)(n.li,{children:"Monoclinic"}),"\n",(0,t.jsx)(n.li,{children:"Orthorhombic"}),"\n",(0,t.jsx)(n.li,{children:"Tetragonal"}),"\n",(0,t.jsx)(n.li,{children:"Trigonal"}),"\n",(0,t.jsx)(n.li,{children:"Hexagonal"}),"\n",(0,t.jsx)(n.li,{children:"Cubic"}),"\n"]}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"str"})," - The name of the crystal system."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_bravais_lattice",children:"get_bravais_lattice"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_bravais_lattice()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Return Bravais lattice in the Pearson notation, where the first\nlowercase letter indicates the crystal system, and the second uppercase\nletter indicates the centring type."}),"\n",(0,t.jsx)(n.p,{children:"Crystal system letters:"}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsx)(n.li,{children:"a = triclinic"}),"\n",(0,t.jsx)(n.li,{children:"m = monoclinic"}),"\n",(0,t.jsx)(n.li,{children:"o = orthorhombic"}),"\n",(0,t.jsx)(n.li,{children:"t = tetragonal"}),"\n",(0,t.jsx)(n.li,{children:"h = hexagonal and trigonal"}),"\n",(0,t.jsx)(n.li,{children:"c = cubic"}),"\n"]}),"\n",(0,t.jsx)(n.p,{children:"Lattice type letters:"}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsx)(n.li,{children:"P = Primitive"}),"\n",(0,t.jsx)(n.li,{children:"S (= A or B or C) = One side/face centred"}),"\n",(0,t.jsx)(n.li,{children:"I = Body centered"}),"\n",(0,t.jsx)(n.li,{children:"R = Rhombohedral centring"}),"\n",(0,t.jsx)(n.li,{children:"F = All faces centred"}),"\n"]}),"\n",(0,t.jsxs)(n.p,{children:[":param"," crystal_system: The crystal system\n",":param"," space_group: The space group number.\n",":type"," crystal_system: str\n",":type"," space_group: int"]}),"\n",(0,t.jsx)(n.p,{children:":return: The Bravais lattice in the Pearson notation.\n:rtype: str"}),"\n",(0,t.jsx)(n.h2,{id:"get_primitive_system",children:"get_primitive_system"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_primitive_system()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Returns a primitive description for this system."}),"\n",(0,t.jsx)(n.p,{children:"This description uses a primitive lattice where positions of the\natoms, and the cell basis vectors are idealized to follow the\nsymmetries that were found with the given precision. This means that\ne.g. the volume, density, angles between basis vectors and basis vector\nlengths may have small deviations from the original system."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"ASE.Atoms"})," - The primitive system."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_conventional_system",children:"get_conventional_system"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_conventional_system()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Used to get the conventional representation of this system."}),"\n",(0,t.jsx)(n.p,{children:"This description uses a conventional lattice where positions of the\natoms, and the cell basis vectors are idealized to follow the\nsymmetries that were found with the given precision. This means that\ne.g. the volume, density, angles between basis vectors and basis vector\nlengths may have small deviations from the original system."}),"\n",(0,t.jsx)(n.h2,{id:"get_rotations",children:"get_rotations"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_rotations()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get the rotational parts of the Seitz matrices that are associated\nwith this space group. Each rotational matrix is accompanied by a\ntranslation with the same index."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"np.ndarray"})," - Rotation matrices."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_translations",children:"get_translations"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_translations()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get the translational parts of the Seitz matrices that are\nassociated with this space group. Each translation is accompanied\nby a rotational matrix with the same index."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"np.ndarray"})," - Translation vectors."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_choice",children:"get_choice"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_choice()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"str"})," - A string specifying the centring, origin and basis vector\nsettings."]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"get_wyckoff_letters_original",children:"get_wyckoff_letters_original"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_wyckoff_letters_original()\n"})}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of str: Wyckoff letters for the atoms in the original system."}),"\n",(0,t.jsx)(n.h2,{id:"get_equivalent_atoms_original",children:"get_equivalent_atoms_original"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_equivalent_atoms_original()\n"})}),"\n",(0,t.jsx)(n.p,{children:"The equivalent atoms are the same as what spglib already outputs, as\nchanges in the wyckoff letters will not afect the equivalence:"}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of int: A list that maps each atom into a symmetry equivalent\nset."}),"\n",(0,t.jsx)(n.h2,{id:"get_wyckoff_letters_conventional",children:"get_wyckoff_letters_conventional"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_wyckoff_letters_conventional()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get the Wyckoff letters of the atoms in the conventional system."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of str: Wyckoff letters."}),"\n",(0,t.jsx)(n.h2,{id:"get_wyckoff_sets_conventional",children:"get_wyckoff_sets_conventional"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_wyckoff_sets_conventional(return_parameters=True)\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get a list of Wyckoff sets for this system. Wyckoff sets combine\ninformation about the atoms and their positions at specific Wyckoff\npositions."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Arguments"}),":"]}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:[(0,t.jsx)(n.code,{children:"return_parameters"})," ",(0,t.jsx)(n.em,{children:"bool"})," - Whether to return the value of possible\nfree Wyckoff parameters. Set to false if they are not needed,\nas their determination can take some time."]}),"\n"]}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsxs)(n.p,{children:["list of WyckoffSets: A list of :class:",(0,t.jsx)(n.code,{children:".WyckoffSet"})," objects for the\nconventional system."]}),"\n",(0,t.jsx)(n.h2,{id:"get_equivalent_atoms_conventional",children:"get_equivalent_atoms_conventional"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_equivalent_atoms_conventional()\n"})}),"\n",(0,t.jsx)(n.p,{children:"List of equivalent atoms in the idealized system."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of int: A list that maps each atom into a symmetry equivalent\nset."}),"\n",(0,t.jsx)(n.h2,{id:"get_wyckoff_letters_primitive",children:"get_wyckoff_letters_primitive"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_wyckoff_letters_primitive()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Get the Wyckoff letters of the atoms in the primitive system."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of str: Wyckoff letters."}),"\n",(0,t.jsx)(n.h2,{id:"get_equivalent_atoms_primitive",children:"get_equivalent_atoms_primitive"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_equivalent_atoms_primitive()\n"})}),"\n",(0,t.jsx)(n.p,{children:"List of equivalent atoms in the primitive system."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:"list of int: A list that maps each atom into a symmetry equivalent\nset."}),"\n",(0,t.jsx)(n.h2,{id:"get_symmetry_dataset",children:"get_symmetry_dataset"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_symmetry_dataset()\n"})}),"\n",(0,t.jsx)(n.p,{children:"Calculates the symmetry dataset with spglib for the given system."}),"\n",(0,t.jsx)(n.h2,{id:"get_symmetry_operations",children:"get_symmetry_operations"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-python",children:"def get_symmetry_operations()\n"})}),"\n",(0,t.jsx)(n.p,{children:"The symmetry operations of the original structure as rotations and\ntranslations."}),"\n",(0,t.jsxs)(n.p,{children:[(0,t.jsx)(n.strong,{children:"Returns"}),":"]}),"\n",(0,t.jsx)(n.p,{children:'Dictionary containing an entry for rotations containing a np.array\nwith 3x3 matrices for each symmetry operation and an entry\n"translations" containing np.array of translations for each\nsymmetry operation.\n3*1 np.ndarray: The shift of the origin as a vector.'})]})}function d(e={}){const{wrapper:n}={...(0,i.a)(),...e.components};return n?(0,t.jsx)(n,{...e,children:(0,t.jsx)(h,{...e})}):h(e)}},1151:(e,n,s)=>{s.d(n,{Z:()=>a,a:()=>l});var t=s(7294);const i={},r=t.createContext(i);function l(e){const n=t.useContext(r);return t.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function a(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(i):e.components||i:l(e.components),t.createElement(r.Provider,{value:n},e.children)}}}]);