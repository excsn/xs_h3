// examples/cell_measures.rs

use xs_h3::{
  cell_area_km2, cell_area_m2, degs_to_rads, get_hexagon_area_avg_km2, get_hexagon_edge_length_avg_km,
  h3_to_string_alloc, lat_lng_to_cell, H3Error, LatLng,
};

fn main() -> Result<(), H3Error> {
  // println!("--- H3 Cell Measures Example ---");

  let point = LatLng {
    lat: degs_to_rads(37.779),
    lng: degs_to_rads(-122.419),
  }; // SF

  for res in (0..=8).step_by(2) {
    // Test a few resolutions
    // println!("\n--- Resolution {} ---", res);

    // Average measures for this resolution
    let avg_area_km2 = get_hexagon_area_avg_km2(res)?;
    let avg_edge_km = get_hexagon_edge_length_avg_km(res)?;
    // println!("Average Hexagon Area at res {}: {:.4} km²", res, avg_area_km2);
    // println!("Average Hexagon Edge Length at res {}: {:.4} km", res, avg_edge_km);

    // Specific cell measures
    let cell = lat_lng_to_cell(&point, res)?;
    // println!("Cell for measures: {}", h3_to_string_alloc(cell));

    let area_km2 = cell_area_km2(cell)?;
    let area_m2 = cell_area_m2(cell)?;
    // println!("Specific Cell Area: {:.4} km² ({:.0} m²)", area_km2, area_m2);

    // Note: Exact edge length requires directed edge H3Indexes.
    // Once that functionality is available, it can be demonstrated here.
    // For now, we only show average edge lengths.
    if xs_h3::is_pentagon(cell) {
      // println!("Note: This cell is a pentagon, exact measures might differ more from averages.");
    }
  }
  Ok(())
}
