//! H3 core library constants.

use std::f64::consts;

use crate::BBox;

// Mathematical constants
/// pi
pub const M_PI: f64 = consts::PI;
/// pi / 2.0
pub const M_PI_2: f64 = consts::FRAC_PI_2;
/// 2.0 * PI
pub const M_2PI: f64 = 2.0 * consts::PI;
/// pi / 180
pub const M_PI_180: f64 = consts::PI / 180.0;
/// 180 / pi
pub const M_180_PI: f64 = 180.0 / consts::PI;

/// Threshold epsilon (general purpose, not tied to a specific unit like degrees or radians)
/// The C library defines this as 0.0000000000000001
pub const EPSILON: f64 = 0.000_000_000_000_000_1;

/// Epsilon for floating point comparisons. ~0.1mm in degrees.
pub const EPSILON_DEG: f64 = 0.000_000_001;
/// Epsilon for floating point comparisons. ~0.1mm in radians.
pub const EPSILON_RAD: f64 = EPSILON_DEG * M_PI_180;

/// sqrt(3) / 2.0, also known as sin(60 degrees)
pub const M_SQRT3_2: f64 = 0.866_025_403_784_438_6; // Equivalent to consts::SQRT_3 / 2.0;

/// Square root of 7.
pub const M_SQRT7: f64 = 2.645_751_311_064_590_6; // sqrt(7.0)
/// Reciprocal of the square root of 7 (1 / sqrt(7)).
pub const M_RSQRT7: f64 = 1.0 / M_SQRT7; // 0.3779644730092272

/// 1 / sin(60 degrees)
pub const M_RSIN60: f64 = 1.0 / M_SQRT3_2;

/// One third
pub const M_ONETHIRD: f64 = 1.0 / 3.0;
/// One seventh
pub const M_ONESEVENTH: f64 = 1.0 / 7.0;

/// Rotation angle between Class II and Class III resolution axes (asin(sqrt(3.0 / 28.0)))
pub const M_AP7_ROT_RADS: f64 = 0.333_473_172_251_832_1;
/// sin(M_AP7_ROT_RADS)
pub const M_SIN_AP7_ROT: f64 = 0.327_326_835_353_988_57;
/// cos(M_AP7_ROT_RADS)
pub const M_COS_AP7_ROT: f64 = 0.944_911_182_523_068_1;

/// Earth radius in kilometers (WGS84 authalic radius)
pub const EARTH_RADIUS_KM: f64 = 6371.007_180_918_475;

/// Scaling factor from hex2d resolution 0 unit length
/// (or distance between adjacent cell center points
/// on the plane) to gnomonic unit length.
pub const RES0_U_GNOMONIC: f64 = 0.381_966_011_250_105; // C: 0.38196601125010500003
/// Inverse of RES0_U_GNOMONIC
pub const INV_RES0_U_GNOMONIC: f64 = 1.0 / RES0_U_GNOMONIC; // C: 2.61803398874989588842

// H3 grid system constants

/// Maximum H3 resolution; H3 has 16 resolutions, numbered 0 through 15.
pub const MAX_H3_RES: i32 = 15;
/// The number of faces on an icosahedron.
pub const NUM_ICOSA_FACES: i32 = 20;
/// The number of H3 base cells.
pub const NUM_BASE_CELLS: i32 = 122;
/// The number of vertices in a hexagon.
pub const NUM_HEX_VERTS: usize = 6;
/// The number of vertices in a pentagon (topologically).
pub const NUM_PENT_VERTS: usize = 5;
/// Maximum number of cell boundary vertices; worst case is a pentagon.
pub const MAX_CELL_BNDRY_VERTS: usize = 10; // 5 original verts + 5 edge crossings

/// The number of pentagons per resolution.
pub const NUM_PENTAGONS: i32 = 12;

// H3 index bit layout constants (as u64 for direct use in bitwise ops)

/// The number of bits in an H3 index.
pub const H3_NUM_BITS: u8 = 64;
/// The bit offset of the max resolution digit in an H3 index (unused).
// pub const H3_MAX_OFFSET: u8 = 63; // This seems to be the high bit itself.

/// The bit offset of the mode in an H3 index.
pub const H3_MODE_OFFSET: u8 = 59;
/// The bit offset of the base cell in an H3 index.
pub const H3_BC_OFFSET: u8 = 45;
/// The bit offset of the resolution in an H3 index.
pub const H3_RES_OFFSET: u8 = 52;
/// The bit offset of the reserved bits in an H3 index.
pub const H3_RESERVED_OFFSET: u8 = 56;
/// The number of bits in a single H3 resolution digit.
pub const H3_PER_DIGIT_OFFSET: u8 = 3;

// Masks for H3 index manipulation
/// 1 in the highest bit, 0's everywhere else.
pub const H3_HIGH_BIT_MASK: u64 = 1u64 << 63;
/// 0 in the highest bit, 1's everywhere else.
pub const H3_HIGH_BIT_MASK_NEGATIVE: u64 = !H3_HIGH_BIT_MASK;
/// 1's in the 4 mode bits, 0's everywhere else.
pub const H3_MODE_MASK: u64 = 0b1111u64 << H3_MODE_OFFSET; // 15
/// 0's in the 4 mode bits, 1's everywhere else.
pub const H3_MODE_MASK_NEGATIVE: u64 = !H3_MODE_MASK;
/// 1's in the 7 base cell bits, 0's everywhere else.
pub const H3_BC_MASK: u64 = 0b111_1111u64 << H3_BC_OFFSET; // 127
/// 0's in the 7 base cell bits, 1's everywhere else.
pub const H3_BC_MASK_NEGATIVE: u64 = !H3_BC_MASK;
/// 1's in the 4 resolution bits, 0's everywhere else.
pub const H3_RES_MASK: u64 = 0b1111u64 << H3_RES_OFFSET; // 15
/// 0's in the 4 resolution bits, 1's everywhere else.
pub const H3_RES_MASK_NEGATIVE: u64 = !H3_RES_MASK;
/// 1's in the 3 reserved bits, 0's everywhere else.
pub const H3_RESERVED_MASK: u64 = 0b111u64 << H3_RESERVED_OFFSET; // 7
/// 0's in the 3 reserved bits, 1's everywhere else.
pub const H3_RESERVED_MASK_NEGATIVE: u64 = !H3_RESERVED_MASK;
/// 1's in the 3 bits of a single H3 digit.
pub const H3_DIGIT_MASK: u64 = 0b111u64; // 7
/// 0's in the 3 bits of a single H3 digit (relative, not used directly as mask often).
// pub const H3_DIGIT_MASK_NEGATIVE: u64 = !H3_DIGIT_MASK; // Be careful with this one.

// H3 index modes
/// Mode for H3 cell indexes.
pub const H3_CELL_MODE: u8 = 1;
/// Mode for H3 directed edge indexes.
pub const H3_DIRECTEDEDGE_MODE: u8 = 2;
/// Mode for H3 undirected edge indexes (currently not a primary H3 type).
pub const H3_EDGE_MODE: u8 = 3;
/// Mode for H3 vertex indexes.
pub const H3_VERTEX_MODE: u8 = 4;

/// H3 index with mode 0, res 0, base cell 0, and 7 for all index digits.
/// Typically used to initialize the creation of an H3 cell index.
pub const H3_INIT: u64 = 35_184_372_088_831;
// More readably:
// H3_INIT = (0b0u64 << 63)                    // High bit
//         | (0b0000u64 << H3_MODE_OFFSET)      // Mode (placeholder, will be set)
//         | (0b000u64 << H3_RESERVED_OFFSET)   // Reserved
//         | (0b0000u64 << H3_RES_OFFSET)       // Resolution (placeholder)
//         | (0b0000000u64 << H3_BC_OFFSET)     // Base cell (placeholder)
//         // Digits (all set to 7, which is 0b111)
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 14))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 13))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 12))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 11))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 10))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 9))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 8))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 7))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 6))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 5))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 4))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 3))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 2))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 1))
//         | (0b111u64 << (H3_PER_DIGIT_OFFSET * 0));
// For H3_INIT, C uses 0x08001fffffffffff which is 35184372088831.
// Let's directly use the value.
// The definition of H3_INIT in C's h3Index.h looks like it's meant to be filled in.
// It has mode 0, res 0, base cell 0, and all 15 digits set to 7 (0b111).
// (0 << 63) | (0 << 59) | (0 << 56) | (0 << 52) | (0 << 45) | (0x7FFF << 0 for digits up to 45 bits)
// If all digits are 7 (0b111), this is 15 * 3 = 45 bits of all 1s.
// (1u64 << 45) - 1
// H3_INIT = 0x0000000000000000 | (0x7fffffffffff); // This would be 45 bits of 1
// Let's use the C literal for now:
// pub const H3_INIT: u64 = 35184372088831;
// The C library's test `setH3Index` implies that H3_INIT already has the digit part set to all 7s.
// Mode 0, Res 0, BC 0, all digits 7.
// Bit layout:
// 1 (reserved) | 4 (mode) | 3 (reserved) | 4 (resolution) | 7 (base cell) | 45 (digits 1-15)
// 0            | 0000     | 000          | 0000           | 0000000       | (15 * 0b111)
// 0x00001FFFFFFFFFFF
// This is `35184372088831`

/// Total number of unique H3 cells at the finest resolution (MAX_H3_RES, which is 15).
/// Formula: 2 + 120 * 7^15
pub const NUM_CELLS_MAX_RES: i64 = 569_707_381_193_162;

/// Maximum cell edge length in radians for each resolution.
/// Computed by taking the max exact edge length for cells
/// at the center of each base cell at that resolution.
#[rustfmt::skip]
pub const MAX_EDGE_LENGTH_RADS: [f64; (MAX_H3_RES + 1) as usize] = [
    0.215_772_062_651_30, // Res 0
    0.083_087_670_684_95, // Res 1
    0.031_489_704_364_39, // Res 2
    0.011_906_628_714_39, // Res 3
    0.004_500_533_309_08, // Res 4
    0.001_701_055_236_19, // Res 5
    0.000_642_939_176_78, // Res 6
    0.000_243_008_206_59, // Res 7
    0.000_091_848_470_87, // Res 8
    0.000_034_715_459_01, // Res 9
    0.000_013_121_210_17, // Res 10
    0.000_004_959_351_29, // Res 11
    0.000_001_874_458_60, // Res 12
    0.000_000_708_478_76, // Res 13
    0.000_000_267_779_80, // Res 14
    0.000_000_101_211_25, // Res 15
];

/// H3 cells that contain the North Pole, by resolution.
#[rustfmt::skip]
pub const NORTH_POLE_CELLS: [u64; (MAX_H3_RES + 1) as usize] = [
    0x8001fffffffffff, 0x81033ffffffffff, 0x820327fffffffff, 0x830326fffffffff,
    0x8403263ffffffff, 0x85032623fffffff, 0x860326237ffffff, 0x870326233ffffff,
    0x880326233bfffff, 0x890326233abffff, 0x8a0326233ab7fff, 0x8b0326233ab0fff,
    0x8c0326233ab03ff, 0x8d0326233ab03bf, 0x8e0326233ab039f, 0x8f0326233ab0399,
];

/// H3 cells that contain the South Pole, by resolution.
#[rustfmt::skip]
pub const SOUTH_POLE_CELLS: [u64; (MAX_H3_RES + 1) as usize] = [
    0x80f3fffffffffff, 0x81f2bffffffffff, 0x82f297fffffffff, 0x83f293fffffffff,
    0x84f2939ffffffff, 0x85f29383fffffff, 0x86f29380fffffff, 0x87f29380effffff,
    0x88f29380e1fffff, 0x89f29380e0fffff, 0x8af29380e0d7fff, 0x8bf29380e0d0fff,
    0x8cf29380e0d0dff, 0x8df29380e0d0cff, 0x8ef29380e0d0cc7, 0x8ff29380e0d0cc4,
];

/// Factor by which to scale the cell bounding box to include all cells.
/// This was determined empirically by finding the smallest factor that
/// passed exhaustive tests.
pub const CELL_SCALE_FACTOR: f64 = 1.1; // From C polyfill.h

/// Factor by which to scale the cell bounding box to include all children.
/// This was determined empirically by finding the smallest factor that
/// passed exhaustive tests.
pub const CHILD_SCALE_FACTOR: f64 = 1.4; // From C polyfill.h

/// Bounding box that defines the valid range for H3 cell centers.
/// Used to quickly check if a LatLng is grossly out of bounds before
/// more expensive H3 operations. {North, South, East, West} in radians.
pub const VALID_RANGE_BBOX: BBox = BBox {
    north: M_PI_2,    // 90 degrees
    south: -M_PI_2,   // -90 degrees
    east: M_PI,       // 180 degrees
    west: -M_PI,      // -180 degrees
};

/// Pre-calculated bounding boxes for all 122 resolution 0 base cells.
/// Values are {north, south, east, west} in radians.
/// Data ported directly from C H3's `polyfill.c`.
#[rustfmt::skip]
pub const RES0_BBOXES: [BBox; NUM_BASE_CELLS as usize] = [
    /*งาม 0 */ BBox { north:  1.5248015836, south:  1.1787242429, east:  2.0562234494, west:  0.4377760900 },
    /*งาม 1 */ BBox { north:  1.5248015836, south:  1.1787242429, east: -0.6066488365, west:  2.5404698032 },
    /*งาม 2 */ BBox { north:  1.5248015836, south:  1.0906938732, east: -2.2499013873, west: -2.8528605332 },
    /*งาม 3 */ BBox { north:  1.4184530255, south:  1.0128514572, east:  0.0056829727, west: -1.1677037961 },
    /*งาม 4 */ BBox { north:  1.2795047789, south:  0.9722665256, east:  0.5555606501, west: -0.1822992482 },
    /*งาม 5 */ BBox { north:  1.3292958660, south:  0.9189892077, east:  2.0562234494, west:  1.0881315428 },
    /*งาม 6 */ BBox { north:  1.3289908609, south:  0.9427181540, east: -2.2987528958, west:  3.0170000806 },
    /*งาม 7 */ BBox { north:  1.2602098388, south:  0.8429122844, east: -0.8997186764, west: -1.7596735929 },
    /*งาม 8 */ BBox { north:  1.2111467388, south:  0.8617060094, east:  1.1912975761, west:  0.4377760900 },
    /*งาม 9 */ BBox { north:  1.2107583144, south:  0.8379533107, east: -1.7202287576, west: -2.4379386170 },
    /*งาม 10*/ BBox { north:  1.1554653095, south:  0.7898245541, east:  2.5365941223, west:  1.8570913345 },
    /*งาม 11*/ BBox { north:  1.1552844509, south:  0.7664142875, east: -3.0673850718, west:  2.5364611027 },
    /*งาม 12*/ BBox { north:  1.1012164356, south:  0.7133009368, east:  0.0964058190, west: -0.5215451449 },
    /*งาม 13*/ BBox { north:  1.0704247279, south:  0.6760394884, east: -0.4798420282, west: -1.1030615958 },
    /*งาม 14*/ BBox { north:  1.0327022877, south:  0.7235635885, east: -2.2499013873, west: -2.7451022089 },
    /*งาม 15*/ BBox { north:  1.0192992465, south:  0.6549123286, east:  0.6303557424, west:  0.0353703010 },
    /*งาม 16*/ BBox { north:  1.0178603759, south:  0.5882763676, east:  1.5319272182, west:  0.9367268251 },
    /*งาม 17*/ BBox { north:  0.9808143416, south:  0.6107606356, east: -2.6710063657, west:  3.0651646303 },
    /*งาม 18*/ BBox { north:  0.9810602322, south:  0.5867983660, east:  2.0282976621, west:  1.5133437497 },
    /*งาม 19*/ BBox { north:  0.9637455181, south:  0.5518649176, east: -1.4297672129, west: -1.9685220251 },
    /*งาม 20*/ BBox { north:  0.8753613623, south:  0.5000895279, east: -1.9243561355, west: -2.4164134319 },
    /*งาม 21*/ BBox { north:  0.8861124347, south:  0.5274296374, east: -0.9578194630, west: -1.4762896628 },
    /*งาม 22*/ BBox { north:  0.8688134327, south:  0.5077056705, east:  1.0323679550, west:  0.5034728403 },
    /*งาม 23*/ BBox { north:  0.8923563821, south:  0.4878126492, east:  2.7643030212, west:  2.2998971670 },
    /*งาม 24*/ BBox { north:  0.8257056928, south:  0.5217310176, east:  2.3092168149, west:  1.9319854185 },
    /*งาม 25*/ BBox { north:  0.8059933046, south:  0.4015081960, east: -3.0641755938, west:  2.7007930081 },
    /*งาม 26*/ BBox { north:  0.8161207973, south:  0.3839680066, east: -0.2161437887, west: -0.7042014970 },
    /*งาม 27*/ BBox { north:  0.7582277987, south:  0.3994355541, east: -2.3405997806, west: -2.8212737380 },
    /*งาม 28*/ BBox { north:  0.7886139100, south:  0.3874201833, east:  0.2311568773, west: -0.2259949106 },
    /*งาม 29*/ BBox { north:  0.7151584037, south:  0.3301247846, east: -0.6484797614, west: -1.0824972810 },
    /*งาม 30*/ BBox { north:  0.7035905107, south:  0.2914867320, east:  1.7144108186, west:  1.2844334838 },
    /*งาม 31*/ BBox { north:  0.6919062957, south:  0.2880831321, east:  0.6486390924, west:  0.1637236928 },
    /*งาม 32*/ BBox { north:  0.6486323568, south:  0.2629042009, east:  2.1031809827, west:  1.6955612255 },
    /*งาม 33*/ BBox { north:  0.6572289230, south:  0.2822265333, east:  1.3091869329, west:  0.8759441627 },
    /*งาม 34*/ BBox { north:  0.6475099776, south:  0.2414986573, east: -1.3027219245, west: -1.6870857014 },
    /*งาม 35*/ BBox { north:  0.6238017405, south:  0.2552208039, east: -2.7242842300, west:  3.1040147326 },
    /*งาม 36*/ BBox { north:  0.6422846044, south:  0.2120675345, east: -1.6763924097, west: -2.1177236674 },
    /*งาม 37*/ BBox { north:  0.5991917539, south:  0.2162046086, east:  2.4859286839, west:  2.0735035389 },
    /*งาม 38*/ BBox { north:  0.5563740687, south:  0.2527655746, east: -0.9988538848, west: -1.3264248933 },
    /*งาม 39*/ BBox { north:  0.5564801333, south:  0.1518740134, east:  2.8703208842, west:  2.4464232048 },
    /*งาม 40*/ BBox { north:  0.5460368800, south:  0.1558909154, east: -2.0678986604, west: -2.4909141961 },
    /*งาม 41*/ BBox { north:  0.5120634778, south:  0.1552202040, east:  0.9544676732, west:  0.5444326211 },
    /*งาม 42*/ BBox { north:  0.4976795156, south:  0.1094489892, east: -0.0433516224, west: -0.4290026815 },
    /*งาม 43*/ BBox { north:  0.4653804551, south:  0.0602996866, east: -0.4124061369, west: -0.8060362378 },
    /*งาม 44*/ BBox { north:  0.4468689109, south:  0.0692685748, east:  0.3205328479, west: -0.0700574888 },
    /*งาม 45*/ BBox { north:  0.4320895823, south:  0.0779644096, east: -3.0623245305, west:  2.8060249999 },
    /*งาม 46*/ BBox { north:  0.4310389261, south:  0.0292743194, east: -2.4158923859, west: -2.8573580993 },
    /*งาม 47*/ BBox { north:  0.3807372758, south: -0.0029701614, east: -0.7703955384, west: -1.1478824872 },
    /*งาม 48*/ BBox { north:  0.3911381671, south: -0.0151876488, east:  1.4913024696, west:  1.1471473174 },
    /*งาม 49*/ BBox { north:  0.3342106317, south:  0.0252661345, east:  1.1514103258, west:  0.8500070626 },
    /*งาม 50*/ BBox { north:  0.3891566980, south: -0.0437135980, east:  1.8804635394, west:  1.4823023138 },
    /*งาม 51*/ BBox { north:  0.3378752085, south: -0.0483509010, east: -1.1227401436, west: -1.4945440882 },
    /*งาม 52*/ BBox { north:  0.3360141896, south: -0.0667506815, east:  2.2379235421, west:  1.8572342301 },
    /*งาม 53*/ BBox { north:  0.3183831810, south: -0.0582195560, east:  0.6605885406, west:  0.2545257294 },
    /*งาม 54*/ BBox { north:  0.3363076150, south: -0.0758954099, east: -1.4795733172, west: -1.8598173569 },
    /*งาม 55*/ BBox { north:  0.2892481735, south: -0.0915063804, east: -1.8356193026, west: -2.2185589736 },
    /*งาม 56*/ BBox { north:  0.2667863228, south: -0.1005808897, east: -2.7680865196, west:  3.1279295327 },
    /*งาม 57*/ BBox { north:  0.2928525414, south: -0.1348316507, east:  2.6140646838, west:  2.2046642291 },
    /*งาม 58*/ BBox { north:  0.2015034281, south: -0.1027985271, east:  0.0688189634, west: -0.2392522941 },
    /*งาม 59*/ BBox { north:  0.2128381330, south: -0.1862683539, east:  2.9380044026, west:  2.5747074766 },
    /*งาม 60*/ BBox { north:  0.1958761421, south: -0.1723703028, east: -2.1694179540, west: -2.5540516588 },
    /*งาม 61*/ BBox { north:  0.1723703033, south: -0.1958761415, east:  0.9721746993, west:  0.5875409945 },
    /*งาม 62*/ BBox { north:  0.1862683544, south: -0.2128381325, east: -0.2035882508, west: -0.5668851768 },
    /*งาม 63*/ BBox { north:  0.1027985275, south: -0.2015034276, east: -3.0727736899, west:  2.9023403595 },
    /*งาม 64*/ BBox { north:  0.1348316512, south: -0.2928525409, east: -0.5275279695, west: -0.9369284242 },
    /*งาม 65*/ BBox { north:  0.1005808902, south: -0.2667863222, east:  0.3735061337, west: -0.0136631208 },
    /*งาม 66*/ BBox { north:  0.0915063809, south: -0.2892481729, east:  1.3059733507, west:  0.9230336798 },
    /*งาม 67*/ BBox { north:  0.0758954106, south: -0.3363076144, east:  1.6620193362, west:  1.2817752964 },
    /*งาม 68*/ BBox { north:  0.0582195565, south: -0.3183831805, east: -2.4810041127, west: -2.8870669240 },
    /*งาม 69*/ BBox { north:  0.0667506820, south: -0.3360141890, east: -0.9036691113, west: -1.2843584232 },
    /*งาม 70*/ BBox { north:  0.0483509015, south: -0.3378752080, east:  2.0188525098, west:  1.6470485652 },
    /*งาม 71*/ BBox { north:  0.0437135985, south: -0.3891566975, east: -1.2611291140, west: -1.6592903395 },
    /*งาม 72*/ BBox { north: -0.0252661340, south: -0.3342106311, east: -1.9901823275, west: -2.2915855907 },
    /*งาม 73*/ BBox { north:  0.0151876493, south: -0.3911381666, east: -1.6502901838, west: -1.9944453360 },
    /*งาม 74*/ BBox { north:  0.0029701618, south: -0.3807372753, east:  2.3711971150, west:  1.9937101662 },
    /*งาม 75*/ BBox { north: -0.0292743189, south: -0.4310389256, east:  0.7257002674, west:  0.2842345541 },
    /*งาม 76*/ BBox { north: -0.0779644091, south: -0.4320895817, east:  0.0792681228, west: -0.3355676534 },
    /*งาม 77*/ BBox { north: -0.0692685743, south: -0.4468689104, east: -2.8210598054, west:  3.0715351648 },
    /*งาม 78*/ BBox { north: -0.0602996861, south: -0.4653804545, east:  2.7291865165, west:  2.3355564155 },
    /*งาม 79*/ BBox { north: -0.1094489886, south: -0.4976795151, east:  3.0982410310, west:  2.7125899718 },
    /*งาม 80*/ BBox { north: -0.1552202035, south: -0.5120634772, east: -2.1871249802, west: -2.5971600322 },
    /*งาม 81*/ BBox { north: -0.1558909148, south: -0.5460368794, east:  1.0736939929, west:  0.6506784573 },
    /*งาม 82*/ BBox { north: -0.1518740130, south: -0.5564801327, east: -0.2712717691, west: -0.6951694486 },
    /*งาม 83*/ BBox { north: -0.2527655741, south: -0.5563740682, east:  2.1427387686, west:  1.8151677600 },
    /*งาม 84*/ BBox { north: -0.2162046081, south: -0.5991917533, east: -0.6556639695, west: -1.0680891144 },
    /*งาม 85*/ BBox { north: -0.2120675340, south: -0.6422846038, east:  1.4652002437, west:  1.0238689859 },
    /*งาม 86*/ BBox { north: -0.2552208034, south: -0.6238017399, east:  0.4173084233, west: -0.0375779209 },
    /*งาม 87*/ BBox { north: -0.2414986568, south: -0.6475099771, east:  1.8388707289, west:  1.4545069520 },
    /*งาม 88*/ BBox { north: -0.2822265329, south: -0.6572289225, east: -1.8324057205, west: -2.2656484906 },
    /*งาม 89*/ BBox { north: -0.2629042004, south: -0.6486323563, east: -1.0384116707, west: -1.4460314278 },
    /*งาม 90*/ BBox { north: -0.2880831316, south: -0.6919062952, east: -2.4929535609, west: -2.9778689605 },
    /*งาม 91*/ BBox { north: -0.2914867316, south: -0.7035905102, east: -1.4271818348, west: -1.8571591695 },
    /*งาม 92*/ BBox { north: -0.3301247841, south: -0.7151584032, east:  2.4931128920, west:  2.0590953724 },
    /*งาม 93*/ BBox { north: -0.3874201828, south: -0.7886139094, east: -2.9104357760, west:  2.9155977430 },
    /*งาม 94*/ BBox { north: -0.3994355536, south: -0.7582277983, east:  0.8009928728, west:  0.3203189154 },
    /*งาม 95*/ BBox { north: -0.3839680061, south: -0.8161207968, east:  2.9254488647, west:  2.4373911564 },
    /*งาม 96*/ BBox { north: -0.4015081955, south: -0.8059933041, east:  0.0774170600, west: -0.4407996455 },
    /*งาม 97*/ BBox { north: -0.5217310172, south: -0.8257056923, east: -0.8323758387, west: -1.2096072351 },
    /*งาม 98*/ BBox { north: -0.4878126487, south: -0.8923563816, east: -0.3772896322, west: -0.8416954863 },
    /*งาม 99*/ BBox { north: -0.5077056699, south: -0.8688134323, east: -2.1092246984, west: -2.6381198131 },
    /*งาม 100*/BBox { north: -0.5274296369, south: -0.8861124342, east:  2.1837731904, west:  1.6653029906 },
    /*งาม 101*/BBox { north: -0.5000895274, south: -0.8753613619, east:  1.2172365179, west:  0.7251792214 },
    /*งาม 102*/BBox { north: -0.5518649171, south: -0.9637455176, east:  1.7118254405, west:  1.1730706283 },
    /*งาม 103*/BBox { north: -0.5867983655, south: -0.9810602317, east: -1.1132949912, west: -1.6282489036 },
    /*งาม 104*/BBox { north: -0.6107606351, south: -0.9808143411, east:  0.4705862876, west: -0.0764280232 },
    /*งาม 105*/BBox { north: -0.5882763671, south: -1.0178603754, east: -1.6096654352, west: -2.2048658282 },
    /*งาม 106*/BBox { north: -0.6549123281, south: -1.0192992459, east: -2.5112369109, west: -3.1062223524 },
    /*งาม 107*/BBox { north: -0.7235635880, south: -1.0327022872, east:  0.8916912664, west:  0.3964904444 },
    /*งาม 108*/BBox { north: -0.6760394879, south: -1.0704247274, east:  2.6617506252, west:  2.0385310576 },
    /*งาม 109*/BBox { north: -0.7133009364, south: -1.1012164351, east: -3.0451868343, west:  2.6200475087 },
    /*งาม 110*/BBox { north: -0.7664142870, south: -1.1552844504, east:  0.0742075819, west: -0.6051315508 },
    /*งาม 111*/BBox { north: -0.7898245536, south: -1.1554653090, east: -0.6049985309, west: -1.2845013188 },
    /*งาม 112*/BBox { north: -0.8379533102, south: -1.2107583139, east:  1.4213638958, west:  0.7036540363 },
    /*งาม 113*/BBox { north: -0.8617060089, south: -1.2111467383, east: -1.9502950772, west: -2.7038165634 },
    /*งาม 114*/BBox { north: -0.8429122839, south: -1.2602098384, east:  2.2418739770, west:  1.3819190605 },
    /*งาม 115*/BBox { north: -0.9427181535, south: -1.3289908604, east:  0.8428397578, west: -0.1245925729 },
    /*งาม 116*/BBox { north: -0.9189892073, south: -1.3292958655, east: -1.0853692039, west: -2.0534611106 },
    /*งาม 117*/BBox { north: -0.9722665251, south: -1.2795047784, east: -2.5860320035, west:  2.9592934054 },
    /*งาม 118*/BBox { north: -1.0128514567, south: -1.4184530251, east: -3.1359096806, west:  1.9738888575 },
    /*งาม 119*/BBox { north: -1.0906938727, south: -1.5248015831, east:  0.2887321209, west: -1.4984857630 },
    /*งาม 120*/BBox { north: -1.1787242424, south: -1.5248015831, east:  2.5349438173, west: -0.6011228503 },
    /*งาม 121*/BBox { north: -1.2030547180, south: -1.5248015831, east: -0.6011228503, west:  2.5349438173 },
];
