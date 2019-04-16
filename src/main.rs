extern crate image;

// use std::str::FromStr;
use image::ColorType;
use image::png::PNGEncoder;
use std::fs::File;
use std::io::Write;

#[allow (dead_code)]

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        writeln!(std::io::stderr(), "Usage: lorenz FILE rho sigma beta").unwrap();
        writeln!(std::io::stderr(), "Example: cargo run lorenz.png 0. 0. 0.").unwrap();
        std::process::exit(1);
    }
    let rho  = 28.; // &args[2].parse().unwrap();  // Parameters for the lorenz attractor
    let sigma = 10.; // &args[3].parse().unwrap();
    let beta = 2.667; // &args[4].parse().unwrap();
    let pic_size = (500 as usize, 500 as usize);
    let number_of_points = 50000 as usize; // Number of points being calculated
    println!("You are using rho = {}, sigma = {} and beta = {} ", rho, sigma, beta);
    // let l_a_points = lorenz(*rho, *sigma, *beta, number_of_points);
    let names = ["pic1.png", "pic2.png", "pic3.png", "pic4.png", "pic5.png", "pic6.png", "pic7.png", "pic8.png", "pic9.png", "pic10.png",
                "pic11.png", "pic12.png", "pic13.png", "pic14.png", "pic15.png", "pic16.png", "pic17.png", "pic18.png", "pic19.png", "pic20.png",
                "pic21.png", "pic22.png", "pic23.png", "pic24.png", "pic25.png", "pic26.png", "pic27.png", "pic28.png", "pic29.png", "pic30.png",
                "pic31.png", "pic32.png", "pic33.png", "pic34.png", "pic35.png", "pic36.png", "pic37.png", "pic38.png", "pic39.png", "pic40.png",
                "pic41.png", "pic42.png", "pic43.png", "pic44.png", "pic45.png", "pic46.png", "pic47.png", "pic48.png", "pic49.png", "pic50.png",
                "pic51.png", "pic52.png", "pic53.png", "pic54.png", "pic55.png", "pic56.png", "pic57.png", "pic58.png", "pic59.png", "pic60.png",
                "pic61.png", "pic62.png", "pic63.png", "pic64.png", "pic65.png", "pic66.png", "pic67.png", "pic68.png", "pic69.png", "pic70.png",
                "pic71.png", "pic72.png", "pic73.png", "pic74.png", "pic75.png", "pic76.png", "pic77.png", "pic78.png", "pic79.png", "pic80.png"];

    // Camera settings
    let mut az = 0.; // Azimut: angle from x-axis arround z-axis -180 to 180
    let mut dec = 0.; // Declination: angle determining heigth above xy-plane -90 to 90
    let dist = 100.; // Distance of camera from (0, 0, 0)
    let offset = (*&args[2].parse().unwrap(), *&args[3].parse().unwrap(), *&args[4].parse().unwrap()); // Center of scene: Point (0,0,0) is translated to

    let l_a_points = lorenz(rho, sigma, beta, number_of_points); // Generate set of points in Lorenz attractor
    for n in 0..6 {  // Rotate camera arround the set
        az += 5.;
        let canvas_points = camera(&l_a_points, az, dec, dist, offset);
        let pixels = coor_to_pixels(canvas_points, pic_size);
        write_image(&names[n], &pixels, pic_size).expect("Error writing PNG file");
    }
    // write_image(&args[1], &pixels, pic_size).expect("Error writing PNG file");
}

fn lorenz(rho: f64, sigma: f64, beta: f64, number_of_points: usize) -> Vec<(f64, f64, f64)>{
    // Generate points in the Lorenz attractor
    println!("{} {} {}", rho, sigma, beta);
    let mut l_a_points = Vec::new();
    let mut x0 = 1.0;
    let mut y0 = 1.0;
    let mut z0 = 1.0;
    let t = 0.001;
    for _ in 0..number_of_points {
        let x = x0 + sigma * (y0 - x0) * t;
        let y = y0 + (x0 * (rho - z0) - y0) * t;
        let z = z0 + (x0 * y0 - beta * z0) * t;
        // println!("x{} y{} z{}", x, y, z);
        l_a_points.push((x as f64, y as f64, -z as f64));

        x0 = x;
        y0 = y;
        z0 = z;
    }
    l_a_points
}

fn camera(l_a_points: &Vec<(f64, f64, f64)>, az: f64, dec: f64, dist: f64, offset: (f64, f64, f64)) -> std::vec::Vec<(f64, f64)> {
    let pi = 3.1415926536;
    let dist_eye_canvas = 0.75;
    let mut canvas_points = Vec::new();
    let eyepoint = (dist*(dec*pi/180.).cos()*(az*pi/180.).cos() + offset.0,
                    dist*(dec*pi/180.).cos()*(az*pi/180.).sin() + offset.1,
                    dist*(dec*pi/180.).sin() + offset.2);
    let to_origo = mult_scalar(-1., eyepoint);
    println!("Eye point at {} {}  {}", eyepoint.0, eyepoint.1, eyepoint.2);
    let canvas_x_axis = norm_vec((-eyepoint.1, eyepoint.0, 0.));
    let canvas_y_axis = crossp(canvas_x_axis, norm_vec(to_origo));
    let canvas_origo = add_vec(add_vec(add_vec(mult_scalar(1. - dist_eye_canvas, eyepoint),
                                       mult_scalar(-0.35*len_vec(eyepoint), canvas_x_axis)),
                                       mult_scalar(-0.35*len_vec(eyepoint), canvas_y_axis)),
                                       offset);
    let value = mult_scalar(1. - dist_eye_canvas, eyepoint);
    println!("canvas dist {}  {} {}", value.0, value.1, value.2);
    let value = mult_scalar(-0.5*len_vec(eyepoint), canvas_x_axis);
    println!("canvas_x_axis {}  {} {}", value.0, value.1, value.2);
    let value = mult_scalar(-0.5*len_vec(eyepoint), canvas_y_axis);
    println!("canvas_y_axis {}  {} {}", value.0, value.1, value.2);
    println!("Canvas origo {} {} {}", canvas_origo.0, canvas_origo.1, canvas_origo.2);

    for n in 0..(l_a_points.len()) {
        let point = l_a_points[n];
        let eye_to_point = ab_vec(eyepoint, point);
        let canvas_x = distance(point, eye_to_point, canvas_origo, canvas_x_axis);
        let canvas_y = distance(point, eye_to_point, canvas_origo, canvas_y_axis);
        // println!("{} {} point {} {} {}", canvas_x, canvas_y, point.0, point.1, point.2);
        canvas_points.push((canvas_x as f64, canvas_y as f64));
    }
    canvas_points
}

fn coor_to_pixels(canvas_points: Vec<(f64, f64)>, pic_size: (usize, usize)) -> Vec<u8> {
    // Project points on the xy-plane
    let mut pixels = vec![255; pic_size.0 * pic_size.1];  // Initialize the array holding the pixels. 255 is black 0 i white
    let x_min = 0.;
    let x_max = 30.;
    let y_min = 10.;
    let y_max = 60.;
    println!("x_min {} x_max {} y_min {} y_max {}", x_min, x_max, y_min, y_max);
    let pic_x = pic_size.0 as f64;
    let pic_y = pic_size.1 as f64;
    for n in 0..(canvas_points.len()) {
        let x = canvas_points[n].0;
        let y = canvas_points[n].1;

        let x_norm = (x - x_min)/(x_max - x_min);  // Move coordinates into first quadrant
        let y_norm = 1. - (y - y_min)/(y_max - y_min);
        let mut position = 1;
        if x_norm < 1. {
            position = (x_norm * pic_x + ((y_norm * pic_y) - ((y_norm * pic_y) % 1.)) * pic_y) as u32; // Translate coordinate to pixel position
        }
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

fn distance(point1: (f64, f64, f64), vec1: (f64, f64, f64), point2: (f64, f64, f64), vec2: (f64, f64, f64)) -> f64 {
    // Return the distance between the line defined by point1 and vec1 and the line defined by point2 and vec2
    let n_vec = crossp(vec1, vec2);
    let distance = dotp(n_vec, ab_vec(point1, point2)).abs()/len_vec(n_vec);
    distance
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

#[test]
fn test_distance() {
    assert_eq!(distance((2., 5., 3.), (7., 2., 3.), (7., 2., 3.), (8., 5., 1.)), 116.0_f64/819.0_f64.sqrt());
}

#[test]
fn test_crossp() {
    assert_eq!(crossp((8., 5., 1.), (7., 2., 3.)), (13., -17., -19.));
}

#[test]
fn test_dotp() {
    assert_eq!(dotp((8., 5., 1.), (7., 2., 3.)), 69.);
}

#[test]
fn test_add_vec() {
    assert_eq!(add_vec((8., 5., 1.), (7., 2., 3.)), (15., 7., 4.));
}

#[test]
fn test_ab_vec() {
    assert_eq!(ab_vec((8., 5., 1.), (7., 2., 3.)),  (-1., -3., 2.));
}

#[test]
fn test_mult_scalar() {
    assert_eq!(mult_scalar(2.0, (1., 1., 1.)), (2., 2., 2.));
    }

#[test]
fn test_len_vec() {
    assert_eq!(len_vec((-1./6.0_f64.sqrt(), 1./6.0_f64.sqrt(), -2./6.0_f64.sqrt())), 1.0);
    assert_eq!(len_vec((-1., 1., -2.)), 6.0_f64.sqrt());
}

#[test]
fn test_norm_vec() {
    assert_eq!(norm_vec((1., 0., 0.)), (1., 0., 0.));
    assert_eq!(norm_vec((-1., 1., -2.)), (-1./6.0_f64.sqrt(), 1./6.0_f64.sqrt(), -2./6.0_f64.sqrt()));
}
