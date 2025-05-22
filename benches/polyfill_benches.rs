use criterion::{criterion_group, criterion_main, Criterion, black_box, BatchSize};
use xs_h3::*; // Import all public items

// --- Fixtures for Polyfill Benchmarks ---

// A simple square GeoLoop (longitude needs care for transmeridian cases)
fn create_simple_square_geoloop(center_lat_deg: f64, center_lng_deg: f64, size_deg: f64) -> GeoLoop {
    let half_size = size_deg / 2.0;
    let verts = vec![
        LatLng { lat: degs_to_rads(center_lat_deg - half_size), lng: degs_to_rads(center_lng_deg - half_size) },
        LatLng { lat: degs_to_rads(center_lat_deg - half_size), lng: degs_to_rads(center_lng_deg + half_size) },
        LatLng { lat: degs_to_rads(center_lat_deg + half_size), lng: degs_to_rads(center_lng_deg + half_size) },
        LatLng { lat: degs_to_rads(center_lat_deg + half_size), lng: degs_to_rads(center_lng_deg - half_size) },
    ];
    GeoLoop {
        num_verts: verts.len(),
        verts,
    }
}

// A simple GeoPolygon with one square outer loop and no holes
fn create_simple_geopolygon() -> GeoPolygon {
    let outer_loop = create_simple_square_geoloop(37.77, -122.41, 0.1); // Approx 0.1 degree square
    GeoPolygon {
        geoloop: outer_loop,
        num_holes: 0,
        holes: Vec::new(),
    }
}

// A GeoPolygon with one outer loop and one smaller inner hole
fn create_donut_geopolygon() -> GeoPolygon {
    let outer_loop = create_simple_square_geoloop(37.77, -122.41, 0.1);
    let hole_loop = create_simple_square_geoloop(37.77, -122.41, 0.05); // Smaller, concentric
    GeoPolygon {
        geoloop: outer_loop,
        num_holes: 1,
        holes: vec![hole_loop],
    }
}


// --- Benchmark Functions ---

fn bench_polygon_to_cells(c: &mut Criterion) {
    let simple_polygon = create_simple_geopolygon();
    let donut_polygon = create_donut_geopolygon();
    
    let mut group = c.benchmark_group("polygon_to_cells");

    for res in [6, 8, 10].iter() { // Test a few resolutions
        group.bench_with_input(format!("simple_poly_res_{}", res), &simple_polygon, |b, poly| {
            let max_size_res = max_polygon_to_cells_size(poly, *res, ContainmentMode::Center as u32);
            if max_size_res.is_err() { return; } // Skip if size estimation fails
            let max_size = max_size_res.unwrap() as usize;
            if max_size == 0 { return; } // Skip if no cells expected
            
            let mut out_cells = vec![H3_NULL; max_size];
            b.iter(|| polygon_to_cells(black_box(poly), black_box(*res), black_box(ContainmentMode::Center as u32), black_box(&mut out_cells)));
        });

        group.bench_with_input(format!("donut_poly_res_{}", res), &donut_polygon, |b, poly| {
            let max_size_res = max_polygon_to_cells_size(poly, *res, ContainmentMode::Center as u32);
             if max_size_res.is_err() { return; }
            let max_size = max_size_res.unwrap() as usize;
            if max_size == 0 { return; }
            
            let mut out_cells = vec![H3_NULL; max_size];
            b.iter(|| polygon_to_cells(black_box(poly), black_box(*res), black_box(ContainmentMode::Center as u32), black_box(&mut out_cells)));
        });
    }
    group.finish();
}

fn bench_cells_to_multi_polygon(c: &mut Criterion) {
    let simple_polygon = create_simple_geopolygon();
    let res = 8; // A moderate resolution for polyfill
    let flags = ContainmentMode::Center as u32;

    // Prepare a set of H3 cells from the simple polygon
    let max_size_res = max_polygon_to_cells_size(&simple_polygon, res, flags);
    if max_size_res.is_err() { 
        // println!("Skipping cells_to_multi_polygon bench: max_size estimation failed.");
        return; 
    }
    let max_size = max_size_res.unwrap() as usize;
    if max_size == 0 {
        // println!("Skipping cells_to_multi_polygon bench: no cells for polyfill.");
        return;
    }
    
    let mut h3_set_vec = vec![H3_NULL; max_size];
    if polygon_to_cells(&simple_polygon, res, flags, &mut h3_set_vec).is_err() {
        // println!("Skipping cells_to_multi_polygon bench: polyfill failed.");
        return;
    }
    // Filter out H3_NULL if polygon_to_cells doesn't fill the whole buffer
    let h3_set: Vec<H3Index> = h3_set_vec.into_iter().filter(|&h| h != H3_NULL).collect();
    if h3_set.is_empty() {
        // println!("Skipping cells_to_multi_polygon bench: polyfill resulted in empty set.");
        return;
    }


    c.bench_function("cells_to_multi_polygon_simple", |b| {
        // Cloning h3_set because cells_to_multi_polygon might eventually take &mut or sort
        b.iter_batched(
            || h3_set.clone(), 
            |data| cells_to_multi_polygon(black_box(&data)),
            BatchSize::SmallInput
        )
    });

    // Benchmark with a donut shape (more complex for cellsToMultiPolygon)
    let donut_polygon = create_donut_geopolygon();
    let max_size_donut_res = max_polygon_to_cells_size(&donut_polygon, res, flags);
     if max_size_donut_res.is_err() { 
        // println!("Skipping cells_to_multi_polygon bench (donut): max_size estimation failed.");
        return; 
    }
    let max_size_donut = max_size_donut_res.unwrap() as usize;
     if max_size_donut == 0 {
        // println!("Skipping cells_to_multi_polygon bench (donut): no cells for polyfill.");
        return;
    }
    let mut h3_set_donut_vec = vec![H3_NULL; max_size_donut];
    if polygon_to_cells(&donut_polygon, res, flags, &mut h3_set_donut_vec).is_err() {
         // println!("Skipping cells_to_multi_polygon bench (donut): polyfill failed.");
        return;
    }
    let h3_set_donut: Vec<H3Index> = h3_set_donut_vec.into_iter().filter(|&h| h != H3_NULL).collect();
    if h3_set_donut.is_empty() {
        // println!("Skipping cells_to_multi_polygon bench (donut): polyfill resulted in empty set.");
        return;
    }

    c.bench_function("cells_to_multi_polygon_donut", |b| {
         b.iter_batched(
            || h3_set_donut.clone(), 
            |data| cells_to_multi_polygon(black_box(&data)),
            BatchSize::SmallInput
        )
    });
}


criterion_group!(polyfill_benches, bench_polygon_to_cells, bench_cells_to_multi_polygon);
criterion_main!(polyfill_benches);