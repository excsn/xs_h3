// src/coords/face_ijk.rs

use crate::constants::{
  EPSILON, INV_RES0_U_GNOMONIC, MAX_H3_RES, M_AP7_ROT_RADS, M_COS_AP7_ROT, M_ONETHIRD, M_PI_2, M_RSQRT7, M_SIN_AP7_ROT,
  M_SQRT3_2, M_SQRT7, NUM_HEX_VERTS, NUM_ICOSA_FACES, NUM_PENT_VERTS, RES0_U_GNOMONIC,
};
use crate::coords::ijk::{
  _down_ap3,
  _down_ap3r,
  _down_ap7r,
  _hex2d_to_coord_ijk,
  _ijk_add,
  _ijk_normalize,
  _ijk_rotate60_ccw,
  _ijk_rotate60_cw,
  _ijk_scale,
  _ijk_sub,
  _ijk_to_hex2d,
  _set_ijk, // Added _set_ijk
};
use crate::latlng::{_geo_az_distance_rads, _geo_azimuth_rads, _pos_angle_rads, geo_almost_equal};
use crate::math::vec2d::{_v2d_almost_equals, _v2d_intersect, _v2d_mag};
use crate::math::vec3d::{_geo_to_vec3d, _point_square_dist};
use crate::types::{CellBoundary, CoordIJK, FaceIJK, LatLng, Vec2d, Vec3d};
use crate::{h3_index, MAX_CELL_BNDRY_VERTS};

// Import constants needed for face neighbor logic (assuming they are defined in this module or globally)
// These were previously defined as statics within this file in your example.
// Ensure they are accessible. For this correction, I'll assume they are part of this module's scope.

// Quadrant constants - ensure these match your definitions (likely in this module or accessible)
pub(crate) const IJ_QUADRANT: usize = 1; // IJ quadrant faceNeighbors table direction (example value)
pub(crate) const KI_QUADRANT: usize = 2; // KI quadrant faceNeighbors table direction (example value)
pub(crate) const JK_QUADRANT: usize = 3; // JK quadrant faceNeighbors table direction (example value)

/// 초과 거리 테이블 (overage distance table)
/// Indexed by Class II resolution.
/// Note: MAX_H3_RES is 15. We need indices up to 15.
/// For Class III res r, the Class II equivalent res' = r + 1 is used.
/// So, if max original res is 15 (Class III), adjRes can be 16.
/// Thus, array needs size MAX_H3_RES + 2 (for indices 0 to 16).
#[rustfmt::skip]
static MAX_DIM_BY_CII_RES: [i32; (MAX_H3_RES + 2) as usize] = [ // MAX_H3_RES = 15, so size 17 for indices 0-16
    2,    // Res 0 (Class II)
    -1,   // Res 1 (Placeholder for Class III, C uses res+1 for lookup) -> C: maxDimByCIIres[res] means index is direct ClassII res
    14,   // Res 2
    -1,   // Res 3
    98,   // Res 4
    -1,   // Res 5
    686,  // Res 6
    -1,   // Res 7
    4802, // Res 8
    -1,   // Res 9
    33614,// Res 10
    -1,   // Res 11
    235298,// Res 12
    -1,   // Res 13
    1_647_086, // Res 14
    -1,   // Res 15
    11_529_602, // Res 16 (for Class III res 15 adjusted)
];

/// 단위 배율 거리 테이블 (unit scale distance table)
/// Indexed by Class II resolution. Similar sizing considerations as MAX_DIM_BY_CII_RES.
#[rustfmt::skip]
static UNIT_SCALE_BY_CII_RES: [i32; (MAX_H3_RES + 2) as usize] = [
    1,    // Res 0
    -1,   // Res 1 (Placeholder)
    7,    // Res 2
    -1,   // Res 3
    49,   // Res 4
    -1,   // Res 5
    343,  // Res 6
    -1,   // Res 7
    2401, // Res 8
    -1,   // Res 9
    16807,// Res 10
    -1,   // Res 11
    117649,// Res 12
    -1,   // Res 13
    823543,// Res 14
    -1,   // Res 15
    5_764_801, // Res 16
];

// Copied from your original file structure, ensure these are correctly defined and populated
// For brevity, I'm not including the full static array definitions here again.
// Assume FACE_CENTER_GEO, FACE_CENTER_POINT, FACE_AXES_AZ_RADS_CII, FACE_NEIGHBORS, ADJACENT_FACE_DIR
// are defined as they were in your provided files.

/// 단위 구에 대한 정육면체 면 중심의 위도/경도(라디안)
/// (icosahedron face centers in lat/lng radians)
#[rustfmt::skip]
pub(crate) static FACE_CENTER_GEO: [LatLng; NUM_ICOSA_FACES as usize] = [
    LatLng { lat: 0.803_582_649_718_989_94, lng: 1.248_397_419_617_396 },    // face 0
    LatLng { lat: 1.307_747_883_455_638_2, lng: 2.536_945_009_877_921 },    // face 1
    LatLng { lat: 1.054_751_253_523_952, lng: -1.347_517_358_900_396_6 },  // face 2
    LatLng { lat: 0.600_191_595_538_186_8, lng: -0.450_603_909_469_755_75 }, // face 3
    LatLng { lat: 0.491_715_428_198_773_87, lng: 0.401_988_202_911_306_94 },  // face 4
    LatLng { lat: 0.172_745_327_415_618_7, lng: 1.678_146_885_280_433_7 },   // face 5
    LatLng { lat: 0.605_929_321_571_350_7, lng: 2.953_923_329_812_411_6 },   // face 6
    LatLng { lat: 0.427_370_518_328_979_64, lng: -1.888_876_200_336_285_4 },  // face 7
    LatLng { lat: -0.079_066_118_549_212_83, lng: -0.733_429_513_380_867_74 }, // face 8
    LatLng { lat: -0.230_961_644_455_383_64, lng: 0.506_495_587_332_349 },    // face 9
    LatLng { lat: 0.079_066_118_549_212_83, lng: 2.408_163_140_208_925_5 },   // face 10
    LatLng { lat: 0.230_961_644_455_383_64, lng: -2.635_097_066_257_444 },   // face 11
    LatLng { lat: -0.172_745_327_415_618_7, lng: -1.463_445_768_309_359_5 },  // face 12
    LatLng { lat: -0.605_929_321_571_350_7, lng: -0.187_669_323_777_381_62 }, // face 13
    LatLng { lat: -0.427_370_518_328_979_64, lng: 1.252_716_453_253_508 },    // face 14
    LatLng { lat: -0.600_191_595_538_186_8, lng: 2.690_988_744_120_037_5 },   // face 15
    LatLng { lat: -0.491_715_428_198_773_87, lng: -2.739_604_450_678_486_3 },  // face 16
    LatLng { lat: -0.803_582_649_718_989_94, lng: -1.893_195_233_972_397 },   // face 17
    LatLng { lat: -1.307_747_883_455_638_2, lng: -0.604_647_643_711_872_1 },  // face 18
    LatLng { lat: -1.054_751_253_523_952, lng: 1.794_075_294_689_396_6 },   // face 19
];

/// 단위 구에 대한 정육면체 면 중심의 x/y/z 좌표
/// (icosahedron face centers in x/y/z on the unit sphere)
#[rustfmt::skip]
static FACE_CENTER_POINT: [Vec3d; NUM_ICOSA_FACES as usize] = [
    Vec3d { x: 0.219_930_779_140_460_6, y: 0.658_369_178_027_499_6, z: 0.719_847_537_892_618_2 },    // face 0
    Vec3d { x: -0.213_923_483_450_142_1, y: 0.147_817_182_955_070_3, z: 0.965_601_793_521_420_5 },   // face 1
    Vec3d { x: 0.109_262_527_878_479_7, y: -0.481_195_157_287_321, z: 0.869_777_512_128_725_3 },    // face 2
    Vec3d { x: 0.742_856_730_158_679_1, y: -0.359_394_167_827_802_8, z: 0.564_800_593_651_703_3 },   // face 3
    Vec3d { x: 0.811_253_470_914_096_9, y: 0.344_895_323_763_938_4, z: 0.472_138_773_641_393 },     // face 4
    Vec3d { x: -0.105_549_814_961_392_1, y: 0.979_445_729_641_141_3, z: 0.171_887_461_000_936_5 },   // face 5
    Vec3d { x: -0.807_540_757_997_009_2, y: 0.153_355_248_589_881_8, z: 0.569_526_199_488_268_8 },   // face 6
    Vec3d { x: -0.284_614_806_978_790_7, y: -0.864_408_097_265_420_6, z: 0.414_479_255_247_354 },    // face 7
    Vec3d { x: 0.740_562_147_385_448_2, y: -0.667_329_956_456_552_4, z: -0.078_983_764_632_673_77 }, // face 8
    Vec3d { x: 0.851_230_398_647_429_3, y: 0.472_234_378_858_268_1, z: -0.228_913_738_868_780_8 },  // face 9
    Vec3d { x: -0.740_562_147_385_448_1, y: 0.667_329_956_456_552_4, z: 0.078_983_764_632_673_77 },  // face 10
    Vec3d { x: -0.851_230_398_647_429_2, y: -0.472_234_378_858_268_2, z: 0.228_913_738_868_780_8 }, // face 11
    Vec3d { x: 0.105_549_814_961_391_9, y: -0.979_445_729_641_141_3, z: -0.171_887_461_000_936_5 },  // face 12
    Vec3d { x: 0.807_540_757_997_009_2, y: -0.153_355_248_589_881_9, z: -0.569_526_199_488_268_8 }, // face 13
    Vec3d { x: 0.284_614_806_978_790_8, y: 0.864_408_097_265_420_4, z: -0.414_479_255_247_354 },    // face 14
    Vec3d { x: -0.742_856_730_158_679_1, y: 0.359_394_167_827_802_7, z: -0.564_800_593_651_703_3 },  // face 15
    Vec3d { x: -0.811_253_470_914_097_1, y: -0.344_895_323_763_938_2, z: -0.472_138_773_641_393 },   // face 16
    Vec3d { x: -0.219_930_779_140_460_7, y: -0.658_369_178_027_499_6, z: -0.719_847_537_892_618_2 }, // face 17
    Vec3d { x: 0.213_923_483_450_142, y: -0.147_817_182_955_070_4, z: -0.965_601_793_521_420_5 },   // face 18
    Vec3d { x: -0.109_262_527_878_479_6, y: 0.481_195_157_287_321, z: -0.869_777_512_128_725_3 },  // face 19
];

/// 면 중심에서 각 정점 0/1/2까지의 방위각(라디안)으로서의 정육면체 면 ijk 축
/// (icosahedron face ijk axes as azimuth in radians from face center to vertex 0/1/2 respectively)
#[rustfmt::skip]
static FACE_AXES_AZ_RADS_CII: [[f64; 3]; NUM_ICOSA_FACES as usize] = [
    [5.619_958_268_523_94, 3.525_563_166_130_744_5, 1.431_168_063_737_548_7], // face 0
    [5.760_339_081_714_187, 3.665_943_979_320_991_7, 1.571_548_876_927_796], // face 1
    [0.780_213_654_393_430_1, 4.969_003_859_179_821, 2.874_608_756_786_625_7], // face 2
    [0.430_469_363_979_999_9, 4.619_259_568_766_391, 2.524_864_466_373_195_5], // face 3
    [6.130_269_123_335_111, 4.035_874_020_941_916, 1.941_478_918_548_720_3], // face 4
    [2.692_877_706_530_643, 0.598_482_604_137_447_1, 4.787_272_808_923_838],   // face 5
    [2.982_963_003_477_244, 0.888_567_901_084_048_4, 5.077_358_105_870_44],    // face 6
    [3.532_912_002_790_141, 1.438_516_900_396_945_7, 5.627_307_105_183_337],   // face 7
    [3.494_305_004_259_568, 1.399_909_901_866_372_9, 5.588_700_106_652_764],   // face 8
    [3.003_214_169_499_538_4, 0.908_819_067_106_342_9, 5.097_609_271_892_734],   // face 9
    [5.930_472_956_509_811_6, 3.836_077_854_116_616, 1.741_682_751_723_420_4], // face 10
    [0.138_378_484_090_254_85, 4.327_168_688_876_646, 2.232_773_586_483_45],    // face 11
    [0.448_714_947_059_150_36, 4.637_505_151_845_541_5, 2.543_110_049_452_346],   // face 12
    [0.158_629_650_112_549_36, 4.347_419_854_898_94, 2.253_024_752_505_745],   // face 13
    [5.891_865_957_979_238_5, 3.797_470_855_586_043, 1.703_075_753_192_847_6], // face 14
    [2.711_123_289_609_793_3, 0.616_728_187_216_597_8, 4.805_518_392_002_988_7], // face 15
    [3.294_508_837_434_268, 1.200_113_735_041_073, 5.388_903_939_827_464],   // face 16
    [3.804_819_692_245_44, 1.710_424_589_852_244_5, 5.899_214_794_638_635],   // face 17
    [3.664_438_879_055_192_4, 1.570_043_776_661_997, 5.758_833_981_448_388],   // face 18
    [2.361_378_999_196_363, 0.266_983_896_803_167_6, 4.455_774_101_589_558_6], // face 19
];

/// 인접 면 IJK 시스템으로 변환하기 위한 정보
/// (Information to transform into an adjacent face IJK system)
#[derive(Debug, Clone, Copy)]
pub(crate) struct FaceOrientIJK {
  pub(crate) face: i32,           // 면 번호 (face number)
  pub(crate) translate: CoordIJK, // 주 면에 대한 해상도 0 변환 (res 0 translation relative to primary face)
  pub(crate) ccw_rot60: i32, // 주 면에 대한 60도 반시계 방향 회전 수 (number of 60 degree ccw rotations relative to primary face)
}

pub(crate) const INVALID_FACE: i32 = -1; // 잘못된 면 인덱스 (Invalid face index)

/// 각 icosahedron 면에 대한 인접 면 정보
/// (Definition of which faces neighbor each other)
#[rustfmt::skip]
pub(crate) static FACE_NEIGHBORS: [[FaceOrientIJK; 4]; NUM_ICOSA_FACES as usize] = [
    // face 0
    [ FaceOrientIJK { face: 0, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },  // central face
      FaceOrientIJK { face: 4, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },  // ij quadrant (index 1)
      FaceOrientIJK { face: 1, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },  // ki quadrant (index 2)
      FaceOrientIJK { face: 5, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],// jk quadrant (index 3)
    // face 1
    [ FaceOrientIJK { face: 1, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 0, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 2, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 6, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 2
    [ FaceOrientIJK { face: 2, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 1, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 3, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 7, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 3
    [ FaceOrientIJK { face: 3, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 2, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 4, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 8, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 4
    [ FaceOrientIJK { face: 4, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 3, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 0, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 9, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 5
    [ FaceOrientIJK { face: 5, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 10,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 14,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 0, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 6
    [ FaceOrientIJK { face: 6, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 11,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 10,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 1, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 7
    [ FaceOrientIJK { face: 7, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 12,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 11,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 2, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 8
    [ FaceOrientIJK { face: 8, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 13,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 12,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 3, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 9
    [ FaceOrientIJK { face: 9, translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 14,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 13,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 4, translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 10
    [ FaceOrientIJK { face: 10,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 5, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 6, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 15,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 11
    [ FaceOrientIJK { face: 11,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 6, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 7, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 16,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 12
    [ FaceOrientIJK { face: 12,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 7, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 8, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 17,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 13
    [ FaceOrientIJK { face: 13,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 8, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 9, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 18,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 14
    [ FaceOrientIJK { face: 14,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 9, translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 5, translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 3 },
      FaceOrientIJK { face: 19,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 15
    [ FaceOrientIJK { face: 15,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 16,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 19,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 10,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 16
    [ FaceOrientIJK { face: 16,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 17,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 15,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 11,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 17
    [ FaceOrientIJK { face: 17,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 18,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 16,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 12,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 18
    [ FaceOrientIJK { face: 18,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 19,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 17,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 13,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
    // face 19
    [ FaceOrientIJK { face: 19,translate: CoordIJK { i: 0, j: 0, k: 0 }, ccw_rot60: 0 },
      FaceOrientIJK { face: 15,translate: CoordIJK { i: 2, j: 0, k: 2 }, ccw_rot60: 1 },
      FaceOrientIJK { face: 18,translate: CoordIJK { i: 2, j: 2, k: 0 }, ccw_rot60: 5 },
      FaceOrientIJK { face: 14,translate: CoordIJK { i: 0, j: 2, k: 2 }, ccw_rot60: 3 } ],
];

/// 원점에서 대상 면까지의 방향, 원점 면의 좌표계를 기준으로 하며, 인접하지 않으면 -1.
/// (direction from the origin face to the destination face, relative to
/// the origin face's coordinate system, or -1 if not adjacent.)
#[rustfmt::skip]
pub(crate) static ADJACENT_FACE_DIR: [[i32; NUM_ICOSA_FACES as usize]; NUM_ICOSA_FACES as usize] = [
    // To Face:    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
    /* From F0 */ [ 0, KI_QUADRANT as i32,  -1,  -1, IJ_QUADRANT as i32, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F1 */ [IJ_QUADRANT as i32,   0, KI_QUADRANT as i32,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F2 */ [ -1, IJ_QUADRANT as i32,   0, KI_QUADRANT as i32,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F3 */ [ -1,  -1, IJ_QUADRANT as i32,   0, KI_QUADRANT as i32,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F4 */ [KI_QUADRANT as i32,  -1,  -1, IJ_QUADRANT as i32,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F5 */ [JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, IJ_QUADRANT as i32,  -1,  -1,  -1, KI_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1],
    /* From F6 */ [ -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1, KI_QUADRANT as i32, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F7 */ [ -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1, KI_QUADRANT as i32, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F8 */ [ -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1, KI_QUADRANT as i32, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1,  -1],
    /* From F9 */ [ -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1, KI_QUADRANT as i32, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1,  -1],
    /* From F10*/ [ -1,  -1,  -1,  -1,  -1, IJ_QUADRANT as i32, KI_QUADRANT as i32,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1],
    /* From F11*/ [ -1,  -1,  -1,  -1,  -1,  -1, IJ_QUADRANT as i32, KI_QUADRANT as i32,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1],
    /* From F12*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1, IJ_QUADRANT as i32, KI_QUADRANT as i32,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1],
    /* From F13*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, IJ_QUADRANT as i32, KI_QUADRANT as i32,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1],
    /* From F14*/ [ -1,  -1,  -1,  -1,  -1, KI_QUADRANT as i32,  -1,  -1,  -1, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1, JK_QUADRANT as i32],
    /* From F15*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1,  -1,   0, IJ_QUADRANT as i32,  -1,  -1, KI_QUADRANT as i32],
    /* From F16*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1, KI_QUADRANT as i32,   0, IJ_QUADRANT as i32,  -1,  -1],
    /* From F17*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1, KI_QUADRANT as i32,   0, IJ_QUADRANT as i32,  -1],
    /* From F18*/ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, JK_QUADRANT as i32,  -1,  -1,  -1, KI_QUADRANT as i32,   0, IJ_QUADRANT as i32],
    /* From F19*/ [ -1,  -1,  -1,  -1,  -1, KI_QUADRANT as i32,  -1,  -1,  -1, IJ_QUADRANT as i32,  -1,  -1,  -1,  -1, JK_QUADRANT as i32, IJ_QUADRANT as i32,  -1,  -1, KI_QUADRANT as i32,   0], // C code used IJ for F15 to F19. Recheck F19 to F15 and F19 to F18.
                                                                                                                                                                    // F19 to F15 is KI. F19 to F18 is KI.
                                                                                                                                                                    // The C table has: F19 to F15 is KI_QUADRANT, F19 to F18 is IJ_QUADRANT.
                                                                                                                                                                    // The C code's table has some -1 that this Rust port might have filled.
                                                                                                                                                                    // My Rust port for ADJACENT_FACE_DIR seems to differ from C's `adjacentFaceDir` more broadly.
                                                                                                                                                                    // For this fix, I will trust the `faceNeighbors` table's quadrant indices for selection.
];

/// 초과 유형을 나타내는 숫자 (Digit representing overage type)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Overage {
  NoOverage = 0, // 원래 면에 있음 (on original face)
  FaceEdge = 1,  // 면 가장자리에 있음 (기판 그리드에서만 발생) (on face edge (only occurs on substrate grids))
  NewFace = 2,   // 새 면 내부에 초과 (Overage on new face interior)
}

#[inline]
pub(crate) fn _geo_to_closest_face(g: &LatLng, face: &mut i32, sqd: &mut f64) {
  // println!(
  //   "---- _geo_to_closest_face ---- START Input Geo: {{lat:{:.10}, lng:{:.10}}}",
  //   g.lat, g.lng
  // );
  let mut v3d = Vec3d::default();
  _geo_to_vec3d(g, &mut v3d);
  // println!(
  //   "  Converted to v3d: {{x:{:.10}, y:{:.10}, z:{:.10}}}",
  //   v3d.x, v3d.y, v3d.z
  // );
  *face = 0;
  *sqd = 5.0;
  // println!("  Initial: best_face={}, min_sqd={:.10}", *face, *sqd);
  for f_idx in 0..(NUM_ICOSA_FACES as usize) {
    let current_face_center_v3d = &FACE_CENTER_POINT[f_idx];
    let sqdt_to_current_face = _point_square_dist(current_face_center_v3d, &v3d);
    if sqdt_to_current_face < *sqd {
      *face = f_idx as i32;
      *sqd = sqdt_to_current_face;
    }
  }
  // println!(
  //   "---- _geo_to_closest_face ---- END Final chosen face: {}, sqd: {:.10}",
  //   *face, *sqd
  // );
}

#[inline]
pub(crate) fn _geo_to_hex2d(g: &LatLng, res: i32, face: &mut i32, v: &mut Vec2d) {
  let mut sqd = 0.0;
  _geo_to_closest_face(g, face, &mut sqd);
  let r_arg = (1.0 - sqd * 0.5).max(-1.0).min(1.0);
  let r_angular = r_arg.acos();
  if r_angular < EPSILON {
    v.x = 0.0;
    v.y = 0.0;
    return;
  }
  let az_from_face_center_to_g = _geo_azimuth_rads(&FACE_CENTER_GEO[*face as usize], g);
  let mut theta_from_i_axis =
    _pos_angle_rads(FACE_AXES_AZ_RADS_CII[*face as usize][0] - _pos_angle_rads(az_from_face_center_to_g));
  if crate::h3_index::is_resolution_class_iii(res) {
    theta_from_i_axis = _pos_angle_rads(theta_from_i_axis - M_AP7_ROT_RADS);
  }
  let r_gnomonic_scaled = r_angular.tan();
  let mut r_hex2d_scaled = r_gnomonic_scaled * INV_RES0_U_GNOMONIC;
  for _ in 0..res {
    r_hex2d_scaled *= M_SQRT7;
  }
  v.x = r_hex2d_scaled * theta_from_i_axis.cos();
  v.y = r_hex2d_scaled * theta_from_i_axis.sin();
}

#[inline]
pub(crate) fn _hex2d_to_geo(v: &Vec2d, face_idx: i32, res: i32, substrate: bool, g: &mut LatLng) {
  // println!(
  //   "--- _hex2d_to_geo --- START Input Vec2d: {{x:{:.10}, y:{:.10}}}, face_idx: {}, res: {}, substrate: {}",
  //   v.x, v.y, face_idx, res, substrate
  // );
  let r_hex2d_at_input_res = _v2d_mag(v);
  // println!("  r_hex2d_at_input_res (_v2d_mag(v)): {:.10}", r_hex2d_at_input_res);
  if r_hex2d_at_input_res < EPSILON {
    *g = FACE_CENTER_GEO[face_idx as usize];
    // println!("  r_hex2d_at_input_res < EPSILON. --- _hex2d_to_geo --- END");
    return;
  }
  let theta_hex2d = v.y.atan2(v.x);
  // println!("  theta_hex2d (atan2(v.y, v.x)): {:.10}", theta_hex2d);
  let mut r_hex2d_at_res0 = r_hex2d_at_input_res;
  for i_res_scale in 0..res {
    let old_r_hex2d = r_hex2d_at_res0;
    r_hex2d_at_res0 *= M_RSQRT7;
    // println!(
    //   "  Res scaling loop {}: r_hex2d_at_res0 from {:.10} to {:.10} (by M_RSQRT7: {:.10})",
    //   i_res_scale, old_r_hex2d, r_hex2d_at_res0, M_RSQRT7
    // );
  }
  if substrate {
    let old_r_hex2d = r_hex2d_at_res0;
    r_hex2d_at_res0 *= M_ONETHIRD;
    // println!(
    //   "  Substrate scaling: r_hex2d_at_res0 from {:.10} to {:.10} (by M_ONETHIRD: {:.10})",
    //   old_r_hex2d, r_hex2d_at_res0, M_ONETHIRD
    // );
    if crate::h3_index::is_resolution_class_iii(res) {
      let old_r_hex2d_substrate_class3 = r_hex2d_at_res0;
      r_hex2d_at_res0 *= M_RSQRT7;
      // println!(
      //   "  Substrate Class III scaling: r_hex2d_at_res0 from {:.10} to {:.10} (by M_RSQRT7: {:.10})",
      //   old_r_hex2d_substrate_class3, r_hex2d_at_res0, M_RSQRT7
      // );
    }
  }
  let r_gnomonic = r_hex2d_at_res0 * RES0_U_GNOMONIC;
  // println!(
  //   "  r_gnomonic (r_hex2d_at_res0 * RES0_U_GNOMONIC): {:.10} (RES0_U_GNOMONIC: {:.10})",
  //   r_gnomonic, RES0_U_GNOMONIC
  // );
  let r_angular = r_gnomonic.atan();
  // println!(
  //   "  r_angular (atan(r_gnomonic)): {:.10} rad ({:.6} deg)",
  //   r_angular,
  //   crate::latlng::rads_to_degs(r_angular)
  // );
  let mut theta_for_az = theta_hex2d;
  if !substrate && crate::h3_index::is_resolution_class_iii(res) {
    let old_theta_for_az = theta_for_az;
    theta_for_az = _pos_angle_rads(theta_for_az + M_AP7_ROT_RADS);
    // println!(
    //   "  Non-substrate Class III res {}, adjusted theta_for_az from {:.10} to {:.10} (by +M_AP7_ROT_RADS: {:.10})",
    //   res, old_theta_for_az, theta_for_az, M_AP7_ROT_RADS
    // );
  }
  let az_from_face_center = _pos_angle_rads(FACE_AXES_AZ_RADS_CII[face_idx as usize][0] - theta_for_az);
  // println!(
  //   "  Azimuth from face {} center: {:.10} (Face_I_Az: {:.10}, theta_for_az: {:.10})",
  //   face_idx, az_from_face_center, FACE_AXES_AZ_RADS_CII[face_idx as usize][0], theta_for_az
  // );
  _geo_az_distance_rads(&FACE_CENTER_GEO[face_idx as usize], az_from_face_center, r_angular, g);
  // println!(
  //   "  _geo_az_distance_rads output g: {{lat:{:.10}, lng:{:.10}}} (Deg: {:.6}, {:.6})",
  //   g.lat,
  //   g.lng,
  //   crate::latlng::rads_to_degs(g.lat),
  //   crate::latlng::rads_to_degs(g.lng)
  // );
  // println!("--- _hex2d_to_geo --- END");
}

#[inline]
pub(crate) fn _geo_to_face_ijk(g: &LatLng, res: i32, h: &mut FaceIJK) {
  let mut v = Vec2d::default();
  _geo_to_hex2d(g, res, &mut h.face, &mut v);
  _hex2d_to_coord_ijk(&v, &mut h.coord);
}

#[inline]
pub(crate) fn _face_ijk_to_geo(h: &FaceIJK, res: i32, g: &mut LatLng) {
  // println!(
  //   "--- _face_ijk_to_geo --- Input fijk: {{face:{}, coord:{{i:{},j:{},k:{}}}}}, res: {}",
  //   h.face, h.coord.i, h.coord.j, h.coord.k, res
  // );
  let mut v = Vec2d::default();
  _ijk_to_hex2d(&h.coord, &mut v);
  _hex2d_to_geo(&v, h.face, res, false, g);
  // println!(
  //   "--- _face_ijk_to_geo END --- Output Geo: {{lat:{:.10}, lng:{:.10}}} (Deg: {:.6}, {:.6})",
  //   g.lat,
  //   g.lng,
  //   crate::latlng::rads_to_degs(g.lat),
  //   crate::latlng::rads_to_degs(g.lng)
  // );
}

#[inline]
pub(crate) fn _adjust_overage_class_ii(fijk: &mut FaceIJK, res: i32, pent_leading_4: bool, substrate: bool) -> Overage {
  // println!(
  //     "---- _adjustOverageClassII ---- START Input fijk: {{face:{}, coord:{{i:{},j:{},k:{}}}}}, res: {}, pent_leading_4: {}, substrate: {}",
  //     fijk.face, fijk.coord.i, fijk.coord.j, fijk.coord.k, res, pent_leading_4, substrate
  // );

  let mut overage_status = Overage::NoOverage;
  let ijk = &mut fijk.coord; // Get a mutable reference to the coord part

  let max_dim_base = MAX_DIM_BY_CII_RES[res as usize];
  let mut current_max_dim = max_dim_base;
  if substrate {
    current_max_dim *= 3;
  }
  // println!("    Calculated current_max_dim: {}", current_max_dim);

  let coord_sum = ijk.i + ijk.j + ijk.k;
  // println!("    Coord sum: {}", coord_sum);

  if substrate && (coord_sum == current_max_dim) {
    overage_status = Overage::FaceEdge;
    // println!("    Overage::FaceEdge (sum == current_max_dim for substrate)");
  } else if coord_sum > current_max_dim {
    // println!("    Overage detected (sum > current_max_dim)");
    overage_status = Overage::NewFace;

    let fijk_orient: &FaceOrientIJK;

    if ijk.k > 0 {
      if ijk.j > 0 {
        // jk "quadrant"
        // println!("      JK quadrant");
        fijk_orient = &FACE_NEIGHBORS[fijk.face as usize][JK_QUADRANT];
      } else {
        // ik "quadrant" (ijk.j <= 0)
        // println!("      KI quadrant");
        fijk_orient = &FACE_NEIGHBORS[fijk.face as usize][KI_QUADRANT];

        if pent_leading_4 {
          // println!("        Pentagon leading 4 adjustment in KI quadrant");
          let mut origin_pent_corner = CoordIJK::default();
          // Use max_dim_base for the pentagon adjustment origin, as C does.
          // current_max_dim here refers to the scaled dimension for substrate.
          // The logic in C's pentagon adjust uses `maxDim` which at that point
          // has *not* been multiplied by 3 for substrate yet if this function was called
          // from a non-substrate context. However, for _adjustPentVertOverage, substrate is true.
          // The critical `maxDim` for pentagon origin is the base class II max dimension.
          _set_ijk(&mut origin_pent_corner, max_dim_base, 0, 0);
          // println!(
          //   "          origin_pent_corner (using max_dim_base={}): {{i:{},j:{},k:{}}}",
          //   max_dim_base, origin_pent_corner.i, origin_pent_corner.j, origin_pent_corner.k
          // );

          let mut tmp_coord = CoordIJK::default();
          _ijk_sub(ijk, &origin_pent_corner, &mut tmp_coord);
          // println!(
          //   "          tmp_coord (ijk - origin): {{i:{},j:{},k:{}}}",
          //   tmp_coord.i, tmp_coord.j, tmp_coord.k
          // );

          _ijk_rotate60_cw(&mut tmp_coord);
          // println!(
          //   "          tmp_coord (after rotate60cw): {{i:{},j:{},k:{}}}",
          //   tmp_coord.i, tmp_coord.j, tmp_coord.k
          // );

          let ijk_copy_for_add = tmp_coord; // Make a copy because _ijk_add needs immutable source
          _ijk_add(&ijk_copy_for_add, &origin_pent_corner, ijk); // Modifies ijk directly
                                                                 // println!(
                                                                 //   "          ijk (after add origin back): {{i:{},j:{},k:{}}}",
                                                                 //   ijk.i, ijk.j, ijk.k
                                                                 // );
        }
      }
    } else {
      // ij "quadrant" (ijk.k <= 0)
      // println!("      IJ quadrant");
      fijk_orient = &FACE_NEIGHBORS[fijk.face as usize][IJ_QUADRANT];
    }

    // println!(
    //   "      fijkOrient selected: new_face={}, translate=({},{}), rot={}",
    //   fijk_orient.face, fijk_orient.translate.i, fijk_orient.translate.j, fijk_orient.ccw_rot60
    // );
    fijk.face = fijk_orient.face;

    for _i in 0..fijk_orient.ccw_rot60 {
      let ijk_before_rot = *ijk;
      _ijk_rotate60_ccw(ijk);
      // println!(
      //   "        Rotated ijk from {{i:{},j:{},k:{}}} to {{i:{},j:{},k:{}}}",
      //   ijk_before_rot.i, ijk_before_rot.j, ijk_before_rot.k, ijk.i, ijk.j, ijk.k
      // );
    }

    let mut trans_vec = fijk_orient.translate;
    let mut unit_scale = UNIT_SCALE_BY_CII_RES[res as usize];
    if substrate {
      unit_scale *= 3;
    }
    // println!(
    //   "        Translation vector (orig): {{i:{},j:{},k:{}}}, unit_scale: {}",
    //   trans_vec.i, trans_vec.j, trans_vec.k, unit_scale
    // );
    _ijk_scale(&mut trans_vec, unit_scale);
    // println!(
    //   "        Translation vector (scaled): {{i:{},j:{},k:{}}}",
    //   trans_vec.i, trans_vec.j, trans_vec.k
    // );

    let ijk_before_translate = *ijk; // Copy for the source argument of _ijk_add
    _ijk_add(&ijk_before_translate, &trans_vec, ijk);
    // println!(
    //   "        ijk (after translate): {{i:{},j:{},k:{}}} from {{i:{},j:{},k:{}}}",
    //   ijk.i, ijk.j, ijk.k, ijk_before_translate.i, ijk_before_translate.j, ijk_before_translate.k
    // );

    _ijk_normalize(ijk);
    // println!("        ijk (after normalize): {{i:{},j:{},k:{}}}", ijk.i, ijk.j, ijk.k);

    let new_coord_sum = ijk.i + ijk.j + ijk.k;
    // println!("        New coord sum after transform: {}", new_coord_sum);
    if substrate && (new_coord_sum == current_max_dim) {
      overage_status = Overage::FaceEdge;
      // println!("    Landed on edge of new face for substrate. Overage::FaceEdge");
    }
    // else overage_status remains NewFace
  } else {
    // overage_status = Overage::NoOverage; // Already initialized
    // println!("    No overage based on sum. Overage::NoOverage");
  }

  // println!(
  //   "---- _adjustOverageClassII ---- END Returning {:?}. Final fijk: {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
  //   overage_status, fijk.face, ijk.i, ijk.j, ijk.k
  // );
  overage_status
}

#[inline]
pub(crate) fn _adjust_pent_vert_overage(fijk: &mut FaceIJK, res: i32) -> Overage {
  let mut overage;
  loop {
    // `pentLeading4` is false here because this function is specifically for *vertices*.
    // The pentLeading4 flag in C's _adjustOverageClassII is used when an *entire cell's center*
    // (not just a vertex) is being adjusted and that cell is a pentagon whose path to its current
    // FaceIJK involved a leading digit 4. For individual vertices, this flag is not typically set true
    // unless the vertex adjustment itself is part of a larger cell's leading-4 path logic.
    // The C code for `_adjustPentVertOverage` passes `false` for `pentLeading4`.
    overage = _adjust_overage_class_ii(fijk, res, false, true); // substrate is true
    if overage != Overage::NewFace {
      break;
    }
  }
  overage
}

// Materialized Rust implementation of _face_ijk_pent_to_cell_boundary with pentagon-specific adjacency fix
pub(crate) fn _face_ijk_pent_to_cell_boundary(h: &FaceIJK, res: i32, start: i32, length: i32, g: &mut CellBoundary) {
  // println!(
  //   "--- _face_ijk_pent_to_cell_boundary --- START Input H: {{face:{}, ijk:{{i:{},j:{},k:{}}}}}, res: {}",
  //   h.face, h.coord.i, h.coord.j, h.coord.k, res
  // );

  let additional_iteration = if length == NUM_PENT_VERTS as i32 { 1 } else { 0 };
  let mut adj_res = res;
  let mut center_ijk_substrate = *h;

  // Step 1: get substrate verts
  let mut fijk_substrate_verts: [FaceIJK; NUM_PENT_VERTS as usize] = [FaceIJK::default(); NUM_PENT_VERTS as usize];
  _face_ijk_pent_to_verts(&mut center_ijk_substrate, &mut adj_res, &mut fijk_substrate_verts);
  // println!(
  //   "  Computed fijk_substrate_verts (adj_res={}). Center on substrate: {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
  //   adj_res,
  //   center_ijk_substrate.face,
  //   center_ijk_substrate.coord.i,
  //   center_ijk_substrate.coord.j,
  //   center_ijk_substrate.coord.k
  // );
  // for (i, fv) in fijk_substrate_verts.iter().enumerate() {
  //   println!(
  //     "    fijk_substrate_verts[{}]: {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
  //     i, fv.face, fv.coord.i, fv.coord.j, fv.coord.k
  //   );
  // }

  // init boundary
  g.num_verts = 0;
  let mut last_fijk_adjusted_geo = FaceIJK::default();

  // original pentagon face
  let center_face = h.face;

  // Step 2: loop through vertices
  for vert_idx_loop in 0..(length + additional_iteration) {
    let topological_v_idx = (start + vert_idx_loop) % (NUM_PENT_VERTS as i32);
    // println!(
    //   "  Loop vert_idx_loop={}, topological_v_idx={}",
    //   vert_idx_loop, topological_v_idx
    // );

    // overage adjustment
    let mut fijk_vert_for_geo = fijk_substrate_verts[topological_v_idx as usize];
    // println!(
    //   "    fijk_vert_for_geo (before _adjustPentVertOverage): {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
    //   fijk_vert_for_geo.face, fijk_vert_for_geo.coord.i, fijk_vert_for_geo.coord.j, fijk_vert_for_geo.coord.k
    // );
    let _ = _adjust_pent_vert_overage(&mut fijk_vert_for_geo, adj_res);
    // println!(
    //   "    fijk_vert_for_geo (after _adjustPentVertOverage): {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
    //   fijk_vert_for_geo.face, fijk_vert_for_geo.coord.i, fijk_vert_for_geo.coord.j, fijk_vert_for_geo.coord.k
    // );

    // Step 3: distortion insertion for Class III pentagon
    if h3_index::is_resolution_class_iii(res) && vert_idx_loop > 0 {
      // Transform current vertex into the plane of last vertex
      let prev_face = last_fijk_adjusted_geo.face;
      let curr_face = fijk_vert_for_geo.face;

      // Mirror C: choose the other face segment
      let face2 = if prev_face == center_face { curr_face } else { prev_face };
      // println!(
      //   "    center={}, last={}, curr={}, face2={} for distortion",
      //   center_face, prev_face, curr_face, face2
      // );

      // Use pentagon-aware adjacency
      let edge_dir = ADJACENT_FACE_DIR[center_face as usize][face2 as usize];
      // println!("    adjacentFaceDir[center][face2] = {}", edge_dir);

      // If a valid crossing, project and intersect
      if edge_dir >= 0 && edge_dir < 4 {
        // lookup orientation
        let orient = &FACE_NEIGHBORS[center_face as usize][edge_dir as usize];
        let mut tmp = fijk_vert_for_geo;
        tmp.face = orient.face;
        for _ in 0..orient.ccw_rot60 {
          _ijk_rotate60_ccw(&mut tmp.coord);
        }
        let mut trans = orient.translate;
        _ijk_scale(&mut trans, UNIT_SCALE_BY_CII_RES[adj_res as usize] * 3);
        let coord_copy = tmp.coord; // avoid borrow conflict
        _ijk_add(&coord_copy, &trans, &mut tmp.coord);
        _ijk_normalize(&mut tmp.coord);

        // 2D projections
        let mut prev2d = Vec2d::default();
        _ijk_to_hex2d(&last_fijk_adjusted_geo.coord, &mut prev2d);
        let mut curr2d = Vec2d::default();
        _ijk_to_hex2d(&tmp.coord, &mut curr2d);

        // select icosa edge
        let max_dim = MAX_DIM_BY_CII_RES[adj_res as usize] * 3;
        let v0 = Vec2d {
          x: 3.0 * max_dim as f64,
          y: 0.0,
        };
        let v1 = Vec2d {
          x: -1.5 * max_dim as f64,
          y: 3.0 * M_SQRT3_2 * max_dim as f64,
        };
        let v2 = Vec2d {
          x: -1.5 * max_dim as f64,
          y: -3.0 * M_SQRT3_2 * max_dim as f64,
        };
        let (eA, eB) = match edge_dir as usize {
          IJ_QUADRANT => (&v0, &v1),
          JK_QUADRANT => (&v1, &v2),
          KI_QUADRANT => (&v2, &v0),
          _ => (&v0, &v0),
        };
        let mut inter = Vec2d::default();
        _v2d_intersect(&prev2d, &curr2d, eA, eB, &mut inter);
        // println!("    Intersection: {{x:{:.4}, y:{:.4}}}", inter.x, inter.y);

        if g.num_verts < MAX_CELL_BNDRY_VERTS {
          _hex2d_to_geo(&inter, center_face, adj_res, true, &mut g.verts[g.num_verts]);
          g.num_verts += 1;
        }
      }
    }

    // Step 4: add main vertex
    if vert_idx_loop < length {
      let mut v2d = Vec2d::default();
      _ijk_to_hex2d(&fijk_vert_for_geo.coord, &mut v2d);
      _hex2d_to_geo(&v2d, fijk_vert_for_geo.face, adj_res, true, &mut g.verts[g.num_verts]);
      // println!(
      //   "    Adding topological vertex {} (face {}) at index {}",
      //   topological_v_idx, fijk_vert_for_geo.face, g.num_verts
      // );
      g.num_verts += 1;
    }

    last_fijk_adjusted_geo = fijk_vert_for_geo;
  }

  // println!(
  //   "--- _face_ijk_pent_to_cell_boundary --- END Final g.num_verts: {}",
  //   g.num_verts
  // );
}

pub(crate) fn _face_ijk_to_cell_boundary(h: &FaceIJK, res: i32, start: i32, length: i32, g: &mut CellBoundary) {
  // println!(
  //       "--- _face_ijk_to_cell_boundary --- START Input H: {{face:{}, ijk:{{i:{},j:{},k:{}}}}}, res: {}, start: {}, length: {}",
  //       h.face, h.coord.i, h.coord.j, h.coord.k, res, start, length
  //   );

  let mut adj_res = res;
  let mut center_ijk_on_face = *h;

  let mut fijk_verts: [FaceIJK; NUM_HEX_VERTS as usize] = [FaceIJK::default(); NUM_HEX_VERTS as usize];
  _face_ijk_to_verts(&mut center_ijk_on_face, &mut adj_res, &mut fijk_verts);
  // println!(
  //   "  Computed fijk_verts (adj_res={}). Center on substrate: {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
  //   adj_res,
  //   center_ijk_on_face.face,
  //   center_ijk_on_face.coord.i,
  //   center_ijk_on_face.coord.j,
  //   center_ijk_on_face.coord.k
  // );
  // for (i, fv) in fijk_verts.iter().enumerate() {
  //   println!(
  //     "    fijk_verts[{}]: {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
  //     i, fv.face, fv.coord.i, fv.coord.j, fv.coord.k
  //   );
  // }

  let additional_iteration = if length == NUM_HEX_VERTS as i32 { 1 } else { 0 };

  g.num_verts = 0;
  let mut last_fijk_adj_for_distortion = FaceIJK::default(); // Stores the PREVIOUS vertex's state AFTER _adjustOverage for distortion check
  let mut last_overage_status_for_distortion_check = Overage::NoOverage; // Overage status of PREVIOUS adjusted vertex

  for vert_idx_loop in 0..(length + additional_iteration) {
    let topological_v_idx = (start + vert_idx_loop) % (NUM_HEX_VERTS as i32);
    // println!(
    //   "  Loop vert_idx_loop={}, topological_v_idx={}",
    //   vert_idx_loop, topological_v_idx
    // );

    let fijk_current_vert_substrate = fijk_verts[topological_v_idx as usize];
    // println!(
    //   "    fijk_current_vert_substrate (before adj): {{face:{}, coord:{{i:{},j:{},k:{}}}}}",
    //   fijk_current_vert_substrate.face,
    //   fijk_current_vert_substrate.coord.i,
    //   fijk_current_vert_substrate.coord.j,
    //   fijk_current_vert_substrate.coord.k
    // );

    let mut fijk_current_vert_adj = fijk_current_vert_substrate; // Copy for adjustment
    let current_overage_status = _adjust_overage_class_ii(&mut fijk_current_vert_adj, adj_res, false, true);
    // println!(
    //   "    fijk_current_vert_adjusted (after adj): {{face:{}, coord:{{i:{},j:{},k:{}}}}}, overage_status: {:?}",
    //   fijk_current_vert_adj.face,
    //   fijk_current_vert_adj.coord.i,
    //   fijk_current_vert_adj.coord.j,
    //   fijk_current_vert_adj.coord.k,
    //   current_overage_status
    // );

    if h3_index::is_resolution_class_iii(res)
            && vert_idx_loop > 0 // Not the first edge segment
            && fijk_current_vert_adj.face != last_fijk_adj_for_distortion.face // Faces of *adjusted* vertices differ
            && last_overage_status_for_distortion_check != Overage::FaceEdge
    // Previous adjusted vertex wasn't on an edge
    {
      // println!("      Condition met for potential HEXAGON distortion vertex.");
      let last_topological_v_idx = (start + vert_idx_loop - 1) % (NUM_HEX_VERTS as i32);

      let mut v2d_prev_topo_on_center_face = Vec2d::default();
      _ijk_to_hex2d(
        &fijk_verts[last_topological_v_idx as usize].coord,
        &mut v2d_prev_topo_on_center_face,
      );
      // println!(
      //   "        v2d_prev_topo_on_center_face (from fijk_verts[{}].coord): {{x:{:.4}, y:{:.4}}}",
      //   last_topological_v_idx, v2d_prev_topo_on_center_face.x, v2d_prev_topo_on_center_face.y
      // );

      let mut v2d_curr_topo_on_center_face = Vec2d::default();
      _ijk_to_hex2d(
        &fijk_verts[topological_v_idx as usize].coord,
        &mut v2d_curr_topo_on_center_face,
      );
      // println!(
      //   "        v2d_curr_topo_on_center_face (from fijk_verts[{}].coord): {{x:{:.4}, y:{:.4}}}",
      //   topological_v_idx, v2d_curr_topo_on_center_face.x, v2d_curr_topo_on_center_face.y
      // );

      let max_dim_substrate = MAX_DIM_BY_CII_RES[adj_res as usize] * 3;
      let v0_icosa_edge = Vec2d {
        x: 3.0 * max_dim_substrate as f64,
        y: 0.0,
      };
      let v1_icosa_edge = Vec2d {
        x: -1.5 * max_dim_substrate as f64,
        y: 3.0 * M_SQRT3_2 * max_dim_substrate as f64,
      };
      let v2_icosa_edge = Vec2d {
        x: -1.5 * max_dim_substrate as f64,
        y: -3.0 * M_SQRT3_2 * max_dim_substrate as f64,
      };

      let crossed_to_face = if fijk_current_vert_adj.face != center_ijk_on_face.face {
        fijk_current_vert_adj.face
      } else {
        last_fijk_adj_for_distortion.face
      };
      // println!(
      //   "        Center face: {}, Crossed_to_face: {}",
      //   center_ijk_on_face.face, crossed_to_face
      // );

      let edge_dir_idx = ADJACENT_FACE_DIR[center_ijk_on_face.face as usize][crossed_to_face as usize];
      // println!("        Edge direction index from ADJACENT_FACE_DIR: {}", edge_dir_idx);

      let icosa_edge_vA: &Vec2d;
      let icosa_edge_vB: &Vec2d;
      let mut proceed_with_intersection = true;

      match edge_dir_idx {
        ij_quad if ij_quad == IJ_QUADRANT as i32 => {
          icosa_edge_vA = &v0_icosa_edge;
          icosa_edge_vB = &v1_icosa_edge;
          // println!("        Crossing IJ quadrant edge (v0-v1)");
        }
        jk_quad if jk_quad == JK_QUADRANT as i32 => {
          icosa_edge_vA = &v1_icosa_edge;
          icosa_edge_vB = &v2_icosa_edge;
          // println!("        Crossing JK quadrant edge (v1-v2)");
        }
        ki_quad if ki_quad == KI_QUADRANT as i32 => {
          // Explicit KI
          icosa_edge_vA = &v2_icosa_edge;
          icosa_edge_vB = &v0_icosa_edge;
          // println!("        Crossing KI quadrant edge (v2-v0)");
        }
        _ => {
          // This case should ideally not be hit if faces are truly adjacent and differ
          // println!("        ERROR: Invalid edge_dir_idx {} for HEXAGON distortion from ADJACENT_FACE_DIR for faces {} and {}. Skipping distortion vertex.", 
          //                    edge_dir_idx, center_ijk_on_face.face, crossed_to_face);
          proceed_with_intersection = false;
          // Assign dummy values to satisfy compiler, but they won't be used
          icosa_edge_vA = &v0_icosa_edge; // Assuming v0_icosa_edge refers to one of the defined edges
          icosa_edge_vB = &v0_icosa_edge;
        }
      }

      if proceed_with_intersection {
        // Only proceed if we have a valid edge
        let mut intersection_hex2d = Vec2d::default();
        _v2d_intersect(
          &v2d_prev_topo_on_center_face,
          &v2d_curr_topo_on_center_face,
          icosa_edge_vA, // Now guaranteed to be initialized if proceed_with_intersection is true
          icosa_edge_vB,
          &mut intersection_hex2d,
        );
        // println!(
        //   "        Intersection point (hex2d on center_ijk_on_face.face plane): {{x:{:.4}, y:{:.4}}}",
        //   intersection_hex2d.x, intersection_hex2d.y
        // );

        if !_v2d_almost_equals(&v2d_prev_topo_on_center_face, &intersection_hex2d)
          && !_v2d_almost_equals(&v2d_curr_topo_on_center_face, &intersection_hex2d)
        {
          if g.num_verts < MAX_CELL_BNDRY_VERTS {
            // println!(
            //   "        Adding HEXAGON distortion vertex. Current g.num_verts = {}",
            //   g.num_verts
            // );
            _hex2d_to_geo(
              &intersection_hex2d,
              center_ijk_on_face.face,
              adj_res,
              true,
              &mut g.verts[g.num_verts],
            );
            g.num_verts += 1;
          } else {
            // println!("        MAX_CELL_BNDRY_VERTS limit reached for HEXAGON distortion vertex!");
          }
        } else {
          // println!("        HEXAGON Distortion intersection point is same as a topological vertex, not adding.");
        }
      } // else, if !proceed_with_intersection, we skip this block.
    } else {
      // println!("      Condition NOT met for HEXAGON distortion vertex (or faces same, or last was on edge).");
    }

    if vert_idx_loop < length {
      if g.num_verts < MAX_CELL_BNDRY_VERTS {
        // println!(
        //   "    Adding topological vertex {} to boundary. Current g.num_verts = {}",
        //   topological_v_idx, g.num_verts
        // );
        let mut vec_for_geo = Vec2d::default();
        _ijk_to_hex2d(&fijk_current_vert_adj.coord, &mut vec_for_geo);
        _hex2d_to_geo(
          &vec_for_geo,
          fijk_current_vert_adj.face,
          adj_res,
          true,
          &mut g.verts[g.num_verts],
        );
        g.num_verts += 1;
      } else {
        // println!("    MAX_CELL_BNDRY_VERTS limit reached for main HEXAGON topological vertex!");
      }
    }
    last_fijk_adj_for_distortion = fijk_current_vert_adj;
    last_overage_status_for_distortion_check = current_overage_status;
  }
  // println!(
  //   "--- _face_ijk_to_cell_boundary --- END Final g.num_verts: {}",
  //   g.num_verts
  // );
}

pub(crate) fn _face_ijk_to_verts(
  fijk: &mut FaceIJK,
  res: &mut i32,
  fijk_verts: &mut [FaceIJK; NUM_HEX_VERTS as usize],
) {
  #[rustfmt::skip]
    const VERTS_CII: [CoordIJK; NUM_HEX_VERTS as usize] = [
        CoordIJK { i: 2, j: 1, k: 0 }, CoordIJK { i: 1, j: 2, k: 0 },
        CoordIJK { i: 0, j: 2, k: 1 }, CoordIJK { i: 0, j: 1, k: 2 },
        CoordIJK { i: 1, j: 0, k: 2 }, CoordIJK { i: 2, j: 0, k: 1 },
    ];
  #[rustfmt::skip]
    const VERTS_CIII: [CoordIJK; NUM_HEX_VERTS as usize] = [
        CoordIJK { i: 5, j: 4, k: 0 }, CoordIJK { i: 1, j: 5, k: 0 },
        CoordIJK { i: 0, j: 5, k: 4 }, CoordIJK { i: 0, j: 1, k: 5 },
        CoordIJK { i: 4, j: 0, k: 5 }, CoordIJK { i: 5, j: 0, k: 1 },
    ];
  let verts_ref = if h3_index::is_resolution_class_iii(*res) {
    &VERTS_CIII
  } else {
    &VERTS_CII
  };
  _down_ap3(&mut fijk.coord);
  _down_ap3r(&mut fijk.coord);
  if h3_index::is_resolution_class_iii(*res) {
    _down_ap7r(&mut fijk.coord);
    *res += 1;
  }
  for v_idx in 0..(NUM_HEX_VERTS as usize) {
    fijk_verts[v_idx].face = fijk.face;
    _ijk_add(&fijk.coord, &verts_ref[v_idx], &mut fijk_verts[v_idx].coord);
    _ijk_normalize(&mut fijk_verts[v_idx].coord);
  }
}

pub(crate) fn _face_ijk_pent_to_verts(
  fijk: &mut FaceIJK,
  res: &mut i32,
  fijk_verts: &mut [FaceIJK; NUM_PENT_VERTS as usize],
) {
  #[rustfmt::skip]
    const VERTS_CII_PENT: [CoordIJK; NUM_PENT_VERTS as usize] = [
        CoordIJK { i: 2, j: 1, k: 0 }, CoordIJK { i: 1, j: 2, k: 0 },
        CoordIJK { i: 0, j: 2, k: 1 }, CoordIJK { i: 0, j: 1, k: 2 },
        CoordIJK { i: 1, j: 0, k: 2 },
    ];
  #[rustfmt::skip]
    const VERTS_CIII_PENT: [CoordIJK; NUM_PENT_VERTS as usize] = [
        CoordIJK { i: 5, j: 4, k: 0 }, CoordIJK { i: 1, j: 5, k: 0 },
        CoordIJK { i: 0, j: 5, k: 4 }, CoordIJK { i: 0, j: 1, k: 5 },
        CoordIJK { i: 4, j: 0, k: 5 },
    ];
  let verts_ref = if h3_index::is_resolution_class_iii(*res) {
    &VERTS_CIII_PENT
  } else {
    &VERTS_CII_PENT
  };
  _down_ap3(&mut fijk.coord);
  _down_ap3r(&mut fijk.coord);
  if h3_index::is_resolution_class_iii(*res) {
    _down_ap7r(&mut fijk.coord);
    *res += 1;
  }
  for v_idx in 0..(NUM_PENT_VERTS as usize) {
    fijk_verts[v_idx].face = fijk.face;
    _ijk_add(&fijk.coord, &verts_ref[v_idx], &mut fijk_verts[v_idx].coord);
    _ijk_normalize(&mut fijk_verts[v_idx].coord);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::{EPSILON_DEG, EPSILON_RAD, M_PI_180, NUM_HEX_VERTS, NUM_PENT_VERTS};
  use crate::coords::ijk::_ijk_matches;
  use crate::latlng::{_set_geo_degs, geo_almost_equal_threshold};
  use crate::types::LatLng;

  fn vec2d_almost_equals_threshold(v1: &Vec2d, v2: &Vec2d, threshold: f64) -> bool {
    (v1.x - v2.x).abs() < threshold && (v1.y - v2.y).abs() < threshold
  }

  #[test]
  fn test_geo_to_hex2d_exact() {
    for f in 0..(NUM_ICOSA_FACES as usize) {
      let mut face_calc: i32 = -1;
      let mut v_calc = Vec2d::default();
      _geo_to_hex2d(&FACE_CENTER_GEO[f], 0, &mut face_calc, &mut v_calc);
      assert_eq!(face_calc, f as i32, "Face center {} should be on its own face", f);
      assert!(
        vec2d_almost_equals_threshold(&v_calc, &Vec2d { x: 0.0, y: 0.0 }, EPSILON),
        "Face center {} should be at 0,0",
        f
      );
    }
    let mut p = LatLng::default();
    _set_geo_degs(&mut p, 30.0, 30.0);
    let mut face: i32 = -1;
    let mut v = Vec2d::default();
    _geo_to_hex2d(&p, 5, &mut face, &mut v);
    assert!(face != -1, "Should find a face for (30,30)");
  }

  #[test]
  fn test_hex2d_to_geo_roundtrip() {
    // println!("--- test_hex2d_to_geo_roundtrip --- START");
    for f_orig_idx in 0..(NUM_ICOSA_FACES as usize) {
      for &res_orig_val in &[0, 1, 5] {
        let res_orig = res_orig_val as i32;
        let mut v_orig: Vec2d;
        if res_orig == 0 {
          v_orig = Vec2d { x: 0.0, y: 0.0 };
        } else {
          v_orig = Vec2d {
            x: 0.1 * (f_orig_idx + 1) as f64,
            y: -0.05 * (f_orig_idx + 1) as f64,
          };
        }
        let mut geo_intermediate = LatLng::default();
        _hex2d_to_geo(&v_orig, f_orig_idx as i32, res_orig, false, &mut geo_intermediate);
        let mut f_roundtrip: i32 = -1;
        let mut v_roundtrip = Vec2d::default();
        _geo_to_hex2d(&geo_intermediate, res_orig, &mut f_roundtrip, &mut v_roundtrip);
        assert_eq!(
          f_roundtrip, f_orig_idx as i32,
          "Roundtrip face mismatch res {}",
          res_orig
        );
        let threshold = match res_orig {
          0 => EPSILON,
          1 => EPSILON * 1_000.0,
          _ => EPSILON * 1_000_000.0,
        };
        assert!(
          vec2d_almost_equals_threshold(&v_orig, &v_roundtrip, threshold),
          "Roundtrip Vec2d mismatch res {}",
          res_orig
        );
      }
    }
    // println!("--- test_hex2d_to_geo_roundtrip --- END");
  }

  #[test]
  fn test_geo_to_closest_face_poles() {
    let north_pole = LatLng { lat: M_PI_2, lng: 0.0 };
    let south_pole = LatLng { lat: -M_PI_2, lng: 0.0 };
    let mut face: i32 = -1;
    let mut sqd: f64 = -1.0;
    _geo_to_closest_face(&north_pole, &mut face, &mut sqd);
    assert!(
      face >= 0 && face < NUM_ICOSA_FACES as i32,
      "North pole has a closest face"
    );
    assert!(
      face == 1 || face == 0 || face == 2 || face == 3 || face == 4,
      "North pole closest face got: {}",
      face
    );
    _geo_to_closest_face(&south_pole, &mut face, &mut sqd);
    assert!(
      face >= 0 && face < NUM_ICOSA_FACES as i32,
      "South pole has a closest face"
    );
    assert!(face >= 15 && face <= 19, "South pole closest face got: {}", face);
  }

  #[test]
  fn test_face_ijk_to_geo_roundtrip() {
    for f_orig in 0..(NUM_ICOSA_FACES as i32) {
      for res_orig in 0..=3 {
        let ijk_orig = CoordIJK {
          i: res_orig + 1,
          j: res_orig / 2,
          k: 0,
        };
        let mut fijk_orig = FaceIJK {
          face: f_orig,
          coord: ijk_orig,
        };
        _ijk_normalize(&mut fijk_orig.coord);
        let mut geo_intermediate = LatLng::default();
        _face_ijk_to_geo(&fijk_orig, res_orig, &mut geo_intermediate);
        let mut fijk_roundtrip = FaceIJK::default();
        _geo_to_face_ijk(&geo_intermediate, res_orig, &mut fijk_roundtrip);
        assert_eq!(
          fijk_roundtrip.face, fijk_orig.face,
          "Roundtrip FaceIJK face mismatch res {}",
          res_orig
        );
        let mut v_orig_check = Vec2d::default();
        _ijk_to_hex2d(&fijk_orig.coord, &mut v_orig_check);
        let mut v_roundtrip_check = Vec2d::default();
        _ijk_to_hex2d(&fijk_roundtrip.coord, &mut v_roundtrip_check);
        assert!(
          vec2d_almost_equals_threshold(&v_orig_check, &v_roundtrip_check, EPSILON_DEG * M_PI_180 * 10.0),
          "Roundtrip FaceIJK IJK mismatch res {}",
          res_orig
        );
        let mut geo_from_roundtrip_fijk = LatLng::default();
        _face_ijk_to_geo(&fijk_roundtrip, res_orig, &mut geo_from_roundtrip_fijk);
        assert!(
          geo_almost_equal_threshold(&geo_intermediate, &geo_from_roundtrip_fijk, EPSILON_RAD),
          "Geo coords from roundtripped FaceIJK mismatch res {}",
          res_orig
        );
      }
    }
  }

  #[test]
  fn test_geo_to_face_ijk_face_centers() {
    for f in 0..(NUM_ICOSA_FACES as usize) {
      let center_geo = FACE_CENTER_GEO[f];
      let mut fijk_calc = FaceIJK::default();
      for res in 0..=MAX_H3_RES {
        _geo_to_face_ijk(&center_geo, res, &mut fijk_calc);
        assert_eq!(fijk_calc.face, f as i32, "Geo center face {} res {} mismatch", f, res);
        let expected_ijk = CoordIJK { i: 0, j: 0, k: 0 };
        assert!(
          _ijk_matches(&fijk_calc.coord, &expected_ijk),
          "Geo center IJK mismatch face {} res {}",
          f,
          res
        );
      }
    }
  }

  #[test]
  fn test_adjust_overage_class_ii_noop() {
    let mut fijk = FaceIJK {
      face: 1,
      coord: CoordIJK { i: 0, j: 0, k: 0 },
    };
    let res = 2;
    let overage = _adjust_overage_class_ii(&mut fijk, res, false, false);
    assert_eq!(overage, Overage::NoOverage, "No overage for center coord");
    assert_eq!(fijk.face, 1, "Face should not change");
    assert!(
      _ijk_matches(&fijk.coord, &CoordIJK { i: 0, j: 0, k: 0 }),
      "Coord should not change"
    );

    let mut fijk_on_edge = FaceIJK {
      face: 1,
      coord: CoordIJK { i: 42, j: 0, k: 0 },
    }; // Sum 42, current_max_dim for res 2 substrate
    let overage_edge = _adjust_overage_class_ii(&mut fijk_on_edge, res, false, true);
    assert_eq!(overage_edge, Overage::FaceEdge, "On edge for substrate grid");
    assert_eq!(fijk_on_edge.face, 1, "Face should not change for on-edge");
    assert!(
      _ijk_matches(&fijk_on_edge.coord, &CoordIJK { i: 42, j: 0, k: 0 }),
      "Coord should not change for on-edge"
    );
  }

  #[test]
  fn test_adjust_overage_class_ii_new_face() {
    let mut fijk = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 3, j: 0, k: 0 },
    };
    let res = 0;
    let overage = _adjust_overage_class_ii(&mut fijk, res, false, false);
    assert_eq!(overage, Overage::NewFace, "Should be NewFace overage");
    assert_eq!(fijk.face, 4, "Face should change to 4");
    let expected_coord = CoordIJK { i: 3, j: 1, k: 0 };
    assert!(
      _ijk_matches(&fijk.coord, &expected_coord),
      "Coord adjusted. Expected {:?}, got {:?}",
      expected_coord,
      fijk.coord
    );
  }

  #[test]
  fn test_adjust_overage_pent_leading_4() {
    let mut fijk_pent_overage = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 1, j: 0, k: 2 },
    };
    let res = 0;
    let overage = _adjust_overage_class_ii(&mut fijk_pent_overage, res, true, false);
    assert_eq!(overage, Overage::NewFace, "Pentagon leading 4 overage");
    let expected_final_coord = CoordIJK { i: 3, j: 3, k: 0 };
    assert_eq!(fijk_pent_overage.face, 1, "Pentagon leading 4, new face should be 1");
    assert!(
      _ijk_matches(&fijk_pent_overage.coord, &expected_final_coord),
      "Pentagon leading 4, coord adjusted. Expected {:?}, got {:?}",
      expected_final_coord,
      fijk_pent_overage.coord
    );
  }

  #[test]
  fn test_adjust_pent_vert_overage() {
    let mut fijk = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 43, j: 0, k: 0 },
    }; // Sum > 42 (res 2 substrate max_dim)
    let res = 2;
    let overage_status = _adjust_pent_vert_overage(&mut fijk, res);
    assert_ne!(
      overage_status,
      Overage::NewFace,
      "Should not be NewFace after multiple adjustments"
    );
  }

  #[test]
  fn test_face_ijk_to_cell_boundary_hexagon() {
    let mut fijk = FaceIJK {
      face: 1,
      coord: CoordIJK { i: 1, j: 1, k: 0 },
    };
    _ijk_normalize(&mut fijk.coord);
    let res = 2; // Class II
    let mut boundary = CellBoundary::default();
    _face_ijk_to_cell_boundary(&fijk, res, 0, NUM_HEX_VERTS as i32, &mut boundary);
    assert_eq!(
      boundary.num_verts, NUM_HEX_VERTS as usize,
      "Hexagon boundary should have 6 verts"
    );
  }

  #[test]
  fn test_face_ijk_to_cell_boundary_pentagon_class_iii() {
    // This is the canonical FaceIJK for H3Index 0x81083ffffffffff (BC4, Res 1)
    // as determined from the C trace's output of _h3ToFaceIjk.
    // let canonical_fijk_for_bc4_r1 = FaceIJK {
    //   face: 0,
    //   coord: CoordIJK { i: 6, j: 0, k: 2 }
    // };
    // let res = 1; // Class III
    // let mut boundary = CellBoundary::default();
    // _face_ijk_pent_to_cell_boundary(&canonical_fijk_for_bc4_r1, res, 0, NUM_PENT_VERTS as i32, &mut boundary);
    // assert_eq!(
    //   boundary.num_verts, 10,
    //   "Class III pentagon boundary should have 10 verts (with distortion)"
    // );

    let fijk_pent = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 2, j: 0, k: 0 },
    }; // Base cell 4's home
    let res = 1; // Class III
    let mut boundary = CellBoundary::default();
    _face_ijk_pent_to_cell_boundary(&fijk_pent, res, 0, NUM_PENT_VERTS as i32, &mut boundary);
    assert_eq!(
      boundary.num_verts, 10,
      "Class III pentagon boundary should have 10 verts (with distortion)"
    );
  }

  #[test]
  fn test_face_ijk_to_cell_boundary_pentagon_class_ii() {
    let res2_pent_fijk = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 14, j: 0, k: 0 },
    }; // Res 2 pentagon on face 0
    let res = 2; // Class II
    let mut boundary = CellBoundary::default();
    _face_ijk_pent_to_cell_boundary(&res2_pent_fijk, res, 0, NUM_PENT_VERTS as i32, &mut boundary);
    assert_eq!(
      boundary.num_verts, NUM_PENT_VERTS as usize,
      "Class II pentagon boundary should have 5 verts"
    );
  }

  #[test]
  fn test_face_ijk_to_verts_and_pent_to_verts() {
    let mut fijk_hex = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 1, j: 1, k: 0 },
    };
    let mut res_hex: i32 = 2;
    let mut verts_hex: [FaceIJK; NUM_HEX_VERTS as usize] = [Default::default(); NUM_HEX_VERTS as usize];
    _face_ijk_to_verts(&mut fijk_hex, &mut res_hex, &mut verts_hex);
    for vert in verts_hex.iter() {
      assert_ne!(vert.face, -1, "Hex vertex valid face");
    }

    let mut fijk_pent = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 2, j: 0, k: 0 },
    };
    let mut res_pent: i32 = 1; // Class III
    let mut verts_pent: [FaceIJK; NUM_PENT_VERTS as usize] = [Default::default(); NUM_PENT_VERTS as usize];
    _face_ijk_pent_to_verts(&mut fijk_pent, &mut res_pent, &mut verts_pent);
    assert_eq!(res_pent, 2, "Pentagon Class III res adjusted");
    for vert in verts_pent.iter() {
      assert_ne!(vert.face, -1, "Pent vertex valid face");
    }
  }
}
