// examples/k_ring_example.rs

use xs_h3::{
  degs_to_rads, grid_disk, grid_distance, h3_to_string_alloc, lat_lng_to_cell, max_grid_disk_size, H3Error, H3Index,
  LatLng, H3_NULL,
};

fn main() -> Result<(), H3Error> {
  // println!("--- K-Ring (Grid Disk) Example ---");

  let origin_lat_deg = 37.779;
  let origin_lng_deg = -122.419;
  let res = 7;
  let k = 2; // k-ring of 2

  let origin_point = LatLng {
    lat: degs_to_rads(origin_lat_deg),
    lng: degs_to_rads(origin_lng_deg),
  };

  let origin_cell = lat_lng_to_cell(&origin_point, res)?;
  // println!(
  //   "Origin Cell: {} (Hex: {:x}) at res {}",
  //   h3_to_string_alloc(origin_cell),
  //   origin_cell.0,
  //   res
  // );
  // println!("Calculating k-ring (gridDisk) with k = {}", k);

  // Determine the maximum size needed for the output array
  let max_size = max_grid_disk_size(k)? as usize;
  let mut k_ring_cells = vec![H3_NULL; max_size];

  // Get the k-ring (disk)
  grid_disk(origin_cell, k, &mut k_ring_cells)?;

  // println!("Cells in k-ring (k={}):", k);
  let mut count = 0;
  for cell_h3 in k_ring_cells.iter() {
    if *cell_h3 != H3_NULL {
      count += 1;
      let distance = grid_distance(origin_cell, *cell_h3)?;
      // println!(
      //   "  Cell: {} (Hex: {:x}), Distance from origin: {}",
      //   h3_to_string_alloc(*cell_h3),
      //   cell_h3.0,
      //   distance
      // );
    }
  }
  // println!("Total cells found in k-ring: {}", count);
  if count > max_size {
    // eprintln!("Error: Found more cells than max_grid_disk_size predicted!");
  }

  Ok(())
}
