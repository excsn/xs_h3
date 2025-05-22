use criterion::{black_box, criterion_group, criterion_main, Criterion};
use xs_h3::*; // Import all public items from your xs_h3 crate

// --- Fixtures for benchmarks ---
fn get_fixed_latlng() -> LatLng {
  LatLng {
    lat: degs_to_rads(37.7749),   // San Francisco latitude
    lng: degs_to_rads(-122.4194), // San Francisco longitude
  }
}

fn get_fixed_h3_index_res5() -> H3Index {
  // H3 index for SF at resolution 5
  // You can get this by running lat_lng_to_cell(&get_fixed_latlng(), 5) once
  // For SF (37.7749, -122.4194) at res 5: 0x85283473fffffff
  H3Index(0x85283473fffffff)
}

fn get_fixed_h3_index_res10() -> H3Index {
  // H3 index for SF at resolution 10
  // For SF (37.7749, -122.4194) at res 10: 0x8a2830828767fff
  H3Index(0x8a2830828767fff)
}

// --- Benchmark Functions ---

fn bench_lat_lng_to_cell(c: &mut Criterion) {
  let latlng = get_fixed_latlng();
  let mut group = c.benchmark_group("lat_lng_to_cell");

  for res in [0, 5, 10, 15].iter() {
    group.bench_with_input(format!("res_{}", res), res, |b, &r| {
      b.iter(|| lat_lng_to_cell(black_box(&latlng), black_box(r)));
    });
  }
  group.finish();
}

fn bench_cell_to_lat_lng(c: &mut Criterion) {
  let cell_res5 = get_fixed_h3_index_res5();
  let cell_res10 = get_fixed_h3_index_res10();

  c.benchmark_group("cell_to_lat_lng")
    .bench_function("res_5", |b| b.iter(|| cell_to_lat_lng(black_box(cell_res5))))
    .bench_function("res_10", |b| b.iter(|| cell_to_lat_lng(black_box(cell_res10))));
}

fn bench_cell_to_boundary(c: &mut Criterion) {
  let cell_res5 = get_fixed_h3_index_res5();
  let cell_res10 = get_fixed_h3_index_res10();

  // Pentagons might be slower due to more complex boundary logic
  let pentagon_res5 = H3Index(0x85080003fffffff); // A known pentagon at res 5 (BC 4 center child)

  c.benchmark_group("cell_to_boundary")
    .bench_function("hex_res_5", |b| b.iter(|| cell_to_boundary(black_box(cell_res5))))
    .bench_function("hex_res_10", |b| b.iter(|| cell_to_boundary(black_box(cell_res10))))
    .bench_function("pent_res_5", |b| b.iter(|| cell_to_boundary(black_box(pentagon_res5))));
}

fn bench_is_valid_cell(c: &mut Criterion) {
  let valid_cell = get_fixed_h3_index_res5();
  let invalid_cell_mode = H3Index(0x05283473fffffff); // Mode 0
  let invalid_cell_digit = H3Index(0x852834727fffffff); // Invalid digit D4

  c.benchmark_group("is_valid_cell")
    .bench_function("valid", |b| b.iter(|| is_valid_cell(black_box(valid_cell))))
    .bench_function("invalid_mode", |b| {
      b.iter(|| is_valid_cell(black_box(invalid_cell_mode)))
    })
    .bench_function("invalid_digit", |b| {
      b.iter(|| is_valid_cell(black_box(invalid_cell_digit)))
    });
}

fn bench_get_resolution(c: &mut Criterion) {
  let cell = get_fixed_h3_index_res10();
  c.bench_function("get_resolution", |b| b.iter(|| get_resolution(black_box(cell))));
}

fn bench_is_pentagon(c: &mut Criterion) {
  let hex_cell = get_fixed_h3_index_res5();
  let pent_cell = H3Index(0x85080003fffffff); // BC 4 based pentagon at res 5

  c.benchmark_group("is_pentagon")
    .bench_function("hexagon", |b| b.iter(|| is_pentagon(black_box(hex_cell))))
    .bench_function("pentagon", |b| b.iter(|| is_pentagon(black_box(pent_cell))));
}

// Register benchmark groups
criterion_group!(
  indexing_benches,
  bench_lat_lng_to_cell,
  bench_cell_to_lat_lng,
  bench_cell_to_boundary,
  bench_is_valid_cell,
  bench_get_resolution,
  bench_is_pentagon
);
criterion_main!(indexing_benches);
