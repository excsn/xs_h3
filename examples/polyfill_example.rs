// examples/polyfill_example.rs

use xs_h3::{
  degs_to_rads, max_polygon_to_cells_size, polygon_to_cells, ContainmentMode, GeoLoop, GeoPolygon, H3Error, H3Index,
  LatLng, H3_NULL,
};

fn main() -> Result<(), H3Error> {
  println!("--- Polyfill Example ---");

  // Define a simple square polygon (e.g., around a park in SF)
  // Note: For GeoJSON, the first and last vertex of a loop must be the same.
  // H3 GeoLoop does not require this; it implicitly closes.
  let outer_loop_coords_deg = vec![
    (37.770, -122.440), // SW
    (37.770, -122.435), // SE
    (37.775, -122.435), // NE
    (37.775, -122.440), // NW
  ];

  let outer_loop_verts_rad: Vec<LatLng> = outer_loop_coords_deg
    .iter()
    .map(|(lat, lng)| LatLng {
      lat: degs_to_rads(*lat),
      lng: degs_to_rads(*lng),
    })
    .collect();

  let geoloop = GeoLoop {
    num_verts: outer_loop_verts_rad.len(),
    verts: outer_loop_verts_rad,
  };

  // For this example, no holes.
  let polygon = GeoPolygon {
    geoloop,
    num_holes: 0,
    holes: Vec::new(),
  };

  let res = 10; // Target H3 resolution for polyfill
  let flags = ContainmentMode::Center as u32; // Polyfill cells whose centers are contained

  println!("Polygon defined with {} vertices.", polygon.geoloop.num_verts);
  println!("Target H3 resolution: {}", res);
  println!("Containment mode: Cell centers");

  // 1. Determine max number of cells needed for the output buffer
  let max_cells = max_polygon_to_cells_size(&polygon, res, flags)?;
  println!("Max H3 cells estimated for polyfill: {}", max_cells);

  if max_cells == 0 {
    println!("No cells estimated, polygon might be too small or empty for this resolution.");
    return Ok(());
  }

  // 2. Allocate memory for the output cells
  let mut polyfill_cells = vec![H3_NULL; max_cells as usize];

  // 3. Perform the polyfill operation
  polygon_to_cells(&polygon, res, flags, &mut polyfill_cells)?;

  // 4. Process the results
  println!(
    "H3 cells covering the polygon ({} actual cells found):",
    polyfill_cells.iter().filter(|&&h| h != H3_NULL).count()
  );
  for (i, cell_h3) in polyfill_cells.iter().enumerate() {
    if *cell_h3 != H3_NULL {
      println!(
        "  Cell {}: {} (Hex: {:x})",
        i,
        xs_h3::h3_to_string_alloc(*cell_h3), // Use crate-level import
        cell_h3.0
      );
    } else {
      // If H3_NULL is found before the end, it means fewer cells were needed than max_cells
      // This is expected as max_polygon_to_cells_size is an overestimate.
      println!("  (Output buffer contains H3_NULL from index {} onwards)", i);
      break;
    }
  }

  // Example of using OverlappingBbox for potentially larger coverage
  let flags_overlap_bbox = ContainmentMode::OverlappingBbox as u32;
  let max_cells_overlap = max_polygon_to_cells_size(&polygon, res, flags_overlap_bbox)?;
  if max_cells_overlap > 0 {
    let mut polyfill_cells_overlap = vec![H3_NULL; max_cells_overlap as usize];
    polygon_to_cells(&polygon, res, flags_overlap_bbox, &mut polyfill_cells_overlap)?;
    println!(
      "\nPolyfill with OverlappingBbox mode found {} cells (max est: {}).",
      polyfill_cells_overlap.iter().filter(|&&h| h != H3_NULL).count(),
      max_cells_overlap
    );
  }

  Ok(())
}
