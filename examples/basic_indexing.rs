use xs_h3::{
  cell_to_boundary, cell_to_lat_lng, degs_to_rads, get_base_cell_number, get_resolution, h3_to_string_alloc,
  is_pentagon, is_valid_cell, lat_lng_to_cell, rads_to_degs, H3Error, H3Index, LatLng,
};

fn main() -> Result<(), H3Error> {
  println!("--- Basic H3 Indexing Example ---");

  // 1. Define a LatLng point (e.g., San Francisco City Hall)
  let lat_deg = 37.779265;
  let lng_deg = -122.419277;
  let point = LatLng {
    lat: degs_to_rads(lat_deg),
    lng: degs_to_rads(lng_deg),
  };
  println!("Original Point: Lat {:.6} deg, Lng {:.6} deg", lat_deg, lng_deg);

  // 2. Find the H3 cell containing this point at resolution 9
  let res = 9;
  let cell: H3Index = lat_lng_to_cell(&point, res)?;
  println!(
    "H3 Cell at res {}: {} (Hex: {:x})",
    res,
    h3_to_string_alloc(cell),
    cell.0
  );

  // 3. Validate the cell
  if is_valid_cell(cell) {
    println!("Cell {} is valid.", h3_to_string_alloc(cell));
  } else {
    println!("Cell {} is NOT valid.", h3_to_string_alloc(cell));
    return Ok(()); // Or handle error
  }

  // 4. Get properties of the cell
  println!("Cell Resolution: {}", get_resolution(cell));
  println!("Cell Base Cell Number: {}", get_base_cell_number(cell));
  println!("Is Pentagon: {}", is_pentagon(cell));

  // 5. Get the center of the cell
  let cell_center: LatLng = cell_to_lat_lng(cell)?;
  println!(
    "Cell Center: Lat {:.6} deg, Lng {:.6} deg",
    rads_to_degs(cell_center.lat),
    rads_to_degs(cell_center.lng)
  );

  // 6. Get the boundary of the cell
  let boundary = cell_to_boundary(cell)?;
  println!("Cell Boundary ({} vertices):", boundary.num_verts);
  for i in 0..boundary.num_verts {
    println!(
      "  Vertex {}: Lat {:.6} deg, Lng {:.6} deg",
      i,
      rads_to_degs(boundary.verts[i].lat),
      rads_to_degs(boundary.verts[i].lng)
    );
  }

  println!("\n--- Hierarchy Example ---");
  // 7. Get parent cell at resolution 5
  let parent_res = 5;
  let parent_cell = xs_h3::cell_to_parent(cell, parent_res)?;
  println!(
    "Parent of {} at res {}: {} (Hex: {:x})",
    h3_to_string_alloc(cell),
    parent_res,
    h3_to_string_alloc(parent_cell),
    parent_cell.0
  );

  // 8. Get the center child of the parent at resolution 7
  let center_child_res = 7;
  let center_child = xs_h3::cell_to_center_child(parent_cell, center_child_res)?;
  println!(
    "Center child of {} at res {}: {} (Hex: {:x})",
    h3_to_string_alloc(parent_cell),
    center_child_res,
    h3_to_string_alloc(center_child),
    center_child.0
  );

  Ok(())
}
