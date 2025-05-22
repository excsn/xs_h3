// examples/cell_traversal.rs

use xs_h3::{
  are_neighbor_cells, degs_to_rads, grid_distance, grid_path_cells, grid_path_cells_size, h3_to_string_alloc,
  lat_lng_to_cell, H3Error, H3Index, LatLng, H3_NULL,
};

fn main() -> Result<(), H3Error> {
  println!("--- H3 Cell Traversal Example ---");

  let res = 7;
  let point1_ll = LatLng {
    lat: degs_to_rads(40.6892),
    lng: degs_to_rads(-74.0445),
  }; // Statue of Liberty
  let point2_ll = LatLng {
    lat: degs_to_rads(40.7060),
    lng: degs_to_rads(-73.9969),
  }; // Near Brooklyn Bridge

  let cell1 = lat_lng_to_cell(&point1_ll, res)?;
  let cell2 = lat_lng_to_cell(&point2_ll, res)?;

  println!(
    "Cell 1 (Statue of Liberty area): {} ({:x})",
    h3_to_string_alloc(cell1),
    cell1.0
  );
  println!(
    "Cell 2 (Brooklyn Bridge area): {} ({:x})",
    h3_to_string_alloc(cell2),
    cell2.0
  );

  // 1. Check if they are neighbors (they likely won't be at res 7 for these points)
  let are_neighbors = are_neighbor_cells(cell1, cell2)?;
  println!("Are Cell 1 and Cell 2 neighbors? {}", are_neighbors);

  // 2. Calculate grid distance
  let distance = grid_distance(cell1, cell2)?;
  println!("Grid distance between Cell 1 and Cell 2: {}", distance);

  // 3. Find the grid path between them
  if distance > 0 && distance < 50 {
    // Limit path display for brevity
    let path_size = grid_path_cells_size(cell1, cell2)? as usize;
    let mut path_cells = vec![H3_NULL; path_size];
    grid_path_cells(cell1, cell2, &mut path_cells)?;

    println!("Path from Cell 1 to Cell 2 ({} cells):", path_size);
    for (i, path_cell) in path_cells.iter().enumerate() {
      if *path_cell != H3_NULL {
        println!("  {}: {} ({:x})", i, h3_to_string_alloc(*path_cell), path_cell.0);
      }
    }
  } else if distance == 0 {
    println!("Cells are the same.");
  } else {
    println!(
      "Cells are too far apart to display path in this example (distance: {}).",
      distance
    );
  }

  // 4. Get direct neighbors of Cell 1 (k-ring of 1, excluding self)
  println!("\nDirect neighbors of Cell 1:");
  let max_k1_size = xs_h3::max_grid_disk_size(1)? as usize;
  let mut k1_disk = vec![H3_NULL; max_k1_size];
  xs_h3::grid_disk(cell1, 1, &mut k1_disk)?;
  for neighbor_h3 in k1_disk.iter() {
    if *neighbor_h3 != H3_NULL && *neighbor_h3 != cell1 {
      println!("  Neighbor: {} ({:x})", h3_to_string_alloc(*neighbor_h3), neighbor_h3.0);
    }
  }
  Ok(())
}