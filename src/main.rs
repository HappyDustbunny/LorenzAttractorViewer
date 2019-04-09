extern crate image;

use std::str::FromStr;
use image::ColorType;
use image::png::PNGEncoder;
use std::fs::File;
use std::io::Write;

#[allow (dead_code)]

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        writeln!(std::io::stderr(), "Usage: lorenz FILE rho sigma beta").unwrap();
        writeln!(std::io::stderr(), "Example: {} lorenz.png 28 10 2.667", args[0]).unwrap();
        std::process::exit(1);
    }
    let rho = &args[2].parse().unwrap();
    let sigma = &args[3].parse().unwrap();
    let beta = &args[4].parse().unwrap();
    let pic_size = (500 as usize, 500 as usize);
    let number_of_points = 50000 as usize; // Number of points being calculated
    println!("You are using rho = {}, sigma = {} and beta = {} ", rho, sigma, beta);
    let l_a_points = lorenz(*rho, *sigma, *beta, number_of_points);

    // Camera settings
    let az = 60.; // Azimut: angle from x-axis arround z-axis -180 to 180
    let dec = 30.; // Declination: angle determining heigth above xy-plane -90 to 90
    let dist = 40.; // Distance of camera from (0, 0, 0)
    let mut l_a_points = Vec::new();
    // let mut pixels = vec![0; pic_size.0 * pic_size.1];

    l_a_points = lorenz(*rho, *sigma, *beta, number_of_points);

    // let pixels = coor_to_pixels(l_a_points, number_of_points, pic_size);
    let pixels = camera(l_a_points, az, dec, dist);

    // write_image(&args[1], &pixels, pic_size).expect("Error writing PNG file");
}

fn lorenz(rho: f64, sigma: f64, beta: f64, number_of_points: usize) -> Vec<(f64, f64, f64)>{
    println!("{} {} {}", rho, sigma, beta);
    let mut l_a_points = Vec::new();
    let mut x0 = 1.0;
    let mut y0 = 1.0;
    let mut z0 = 1.0;
    let t = 0.01;
    for _ in 0..number_of_points {
        let x = x0 + sigma * (y0 - x0) * t;
        let y = y0 + (x0 * (rho - z0) - y0) * t;
        let z = z0 + (x0 * y0 - beta * z0) * t;
        // println!("x{} y{} z{}", x, y, z);
        l_a_points.push((x as f64, y as f64, z as f64));

        x0 = x;
        y0 = y;
        z0 = z;
    }
    l_a_points
}

fn camera(l_a_points: Vec<(f64, f64, f64)>, az: f64, dec: f64, dist: f64) -> std::vec::Vec<(f64, f64, f64)> {
    let pi = 3.1415926536;
    let dist_eye_canvas = 0.33;
    let eyepoint = (dist*(dec*pi/180.).cos()*(az*pi/180.).cos(),
                    dist*(dec*pi/180.).cos()*(az*pi/180.).sin(),
                    dist*(dec*pi/180.).sin());
    let to_origo = mult_scalar(-1., eyepoint);
    println!("Eye point at {} {}  {}", eyepoint.0, eyepoint.1, eyepoint.2);
    let canvas_x_axis = norm_vec((-eyepoint.1, eyepoint.0, 0.));
    let canvas_y_axis = crossp(canvas_x, norm_vec(to_origo));
    let canvas_origo = add_vec(add_vec(mult_scalar((1. - dist_eye_canvas), eyepoint),
                                       mult_scalar(-0.5*len_vec(eyepoint), canvas_x)),
                                       mult_scalar(-0.5*len_vec(eyepoint), canvas_y));

    for n in 0..(l_a_points.len()) {
        point = l_a_points[n];
        eye_to_point = ab_vec(eyepoint, point);
        let canvas_x = distace(point, eye_to_point, canvas_origo, canvas_x_axis);
        let canvas_y = distace(point, eye_to_point, canvas_origo, canvas_y_axis);

    }
    println!("Canvas origo {} {} {}", canvas_origo.0, canvas_origo.1, canvas_origo.2);
    l_a_points
}

fn distace(point1: (f64, f64, f64), vec1: (f64, f64, f64), point2: (f64, f64, f64), vec2: (f64, f64, f64)) -> (f64, f64, f64) {
    // Return the distance between the line defined by point1 and vec1 and the line defined by point2 and vec2
    let n_vec = crossp(vec1, vec2);
    let distace = dotp(n_vec, ab_vec(point1, point2))/len_vec(n_vec);
    distace
}


fn crossp(vec1: (f64, f64, f64), vec2: (f64, f64, f64)) -> (f64, f64, f64) {
    let res_vec = (vec1.1*vec2.2 - vec2.1*vec1.2, vec1.2*vec2.0 - vec2.2*vec1.0, vec1.0*vec2.1 - vec2.0*vec1.1);
    res_vec
}

fn dotp(vec1: (f64, f64, f64), vec2: (f64, f64, f64)) -> f64 {
    let dotp = vec1.0*vec2.0 + vec1.1*vec2.1 + vec1.2*vec2.2;
    dotp
}

fn add_vec(vec1: (f64, f64, f64), vec2: (f64, f64, f64)) -> (f64, f64, f64) {
    let res_vec = (vec1.0 + vec2.0, vec1.1 + vec2.1, vec1.2 + vec2.2);
    res_vec
}

fn ab_vec(vec_a: (f64, f64, f64), vec_b: (f64, f64, f64)) -> (f64, f64, f64) {
    // Vector from point A to point B
    let vec_res = (vec_b.0 - vec_a.0, vec_b.1 - vec_a.1, vec_b.2 - vec_a.2);
    vec_res
}

fn mult_scalar(scalar: f64, vec: (f64, f64, f64)) -> (f64, f64, f64) {
    // Multiplies vector vec with a scalar
    let vec_res = (scalar*vec.0, scalar*vec.1, scalar*vec.2);
    vec_res
}

fn len_vec(vec: (f64, f64, f64)) -> f64 {
    // Return length of a vector
    let len_vec = (vec.0*vec.0 + vec.1*vec.1 + vec.2*vec.2).sqrt();
    len_vec
}

fn norm_vec(vec: (f64, f64, f64)) -> (f64, f64, f64) {
    // Return vector of length 1
    let vec_res = mult_scalar(1./len_vec(vec), vec);
    vec_res
}

fn coor_to_pixels(l_a_points: Vec<(f64, f64, f64)>, number_of_points: usize, pic_size: (usize, usize)) -> Vec<u8> {
    // Project points on the xy-plane
    let mut pixels = vec![255; pic_size.0 * pic_size.1];  // Initialize the array holding the pixels. 255 is black 0 i white
    let x_min = -20.;
    let x_max = 20.;
    let y_min = -25.;
    let y_max = 25.;
    let pic_x = pic_size.0 as f64;
    let pic_y = pic_size.1 as f64;
    for n in 0..(number_of_points) {
        let x = l_a_points[n].0;
        let y = l_a_points[n].1;
        let z = l_a_points[n].2;
        let x_norm = (x - x_min)/(x_max - x_min);  // Move coordinates into first quadrant
        let y_norm = 1. - (y - y_min)/(y_max - y_min);
        let position = (x_norm * pic_x + ((y_norm * pic_y) - ((y_norm * pic_y) % 1.)) * pic_y) as u32; // Translate coordinate to pixel position
        // let position = (x/10. * 500. + (((1. - y/10.) * 500.) - (((1. - y/10.) * 500.) % 1.)) * 500.) as u32;
        // println!("x {} y {} z {} pos {}", x, y, z, position);
        if position > 0 && position < 250000 {
            pixels[position as usize] = 0;
        }
    }
    pixels
}

fn write_image(filename: &str, pixels: &[u8], bounds: (usize, usize)) -> Result<(), std::io::Error> {
    let output = File::create(filename)?;
    let encoder = PNGEncoder::new(output);
    encoder.encode(&pixels, bounds.0 as u32, bounds.1 as u32, ColorType::Gray(8))?;
    Ok(())
}
