use kiss3d::nalgebra::{Point3, Translation3};

use kiss3d::light::Light;
use kiss3d::window::Window;

use itertools::izip;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::{thread, time};

fn draw_box(color: &Point3<f32>, x: f32, y: f32, z: f32, window: &mut Window) {
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(x, 0., 0.), color);
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(0., y, 0.), color);
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(0., 0., z), color);
    window.draw_line(&Point3::new(x, 0., 0.), &Point3::new(x, y, 0.), color);
    window.draw_line(&Point3::new(x, 0., 0.), &Point3::new(x, 0., z), color);
    window.draw_line(&Point3::new(x, y, 0.), &Point3::new(0., y, 0.), color);
    window.draw_line(&Point3::new(x, y, 0.), &Point3::new(x, y, z), color);
    window.draw_line(&Point3::new(0., y, 0.), &Point3::new(0., y, z), color);
    window.draw_line(&Point3::new(0., 0., z), &Point3::new(x, 0., z), color);
    window.draw_line(&Point3::new(0., 0., z), &Point3::new(0., y, z), color);
    window.draw_line(&Point3::new(x, 0., z), &Point3::new(x, y, z), color);
    window.draw_line(&Point3::new(x, y, z), &Point3::new(0., y, z), color);
}

fn main() {
    let mut generator = StdRng::seed_from_u64(12);
    let time_step = 1.;
    let max_x: f64 = 60.;
    let max_y: f64 = 60.;
    let max_z: f64 = 60.;
    let max_velocity: f64 = 2.;
    let draw_frequency = time::Duration::from_millis(15);
    let num_spheres = 50;
    let mut sphere_positions_x: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_x)
        .collect();
    let mut sphere_positions_y: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_y)
        .collect();
    let mut sphere_positions_z: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_z)
        .collect();
    let mut sphere_velocities_x: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_velocity)
        .collect();
    let mut sphere_velocities_y: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_velocity)
        .collect();
    let mut sphere_velocities_z: Vec<f64> = (0..num_spheres)
        .map(|_| generator.gen::<f64>() * max_velocity)
        .collect();

    let mut window = Window::new("Kiss3d: cube");
    let mut spheres: Vec<_> = izip!(
        &sphere_positions_x,
        &sphere_positions_y,
        &sphere_positions_z,
    )
    .map(|(x, y, z)| {
        let mut sphere = window.add_sphere(1.0);
        sphere.set_local_translation(Translation3::new(*x as f32, *y as f32, *z as f32));
        sphere.set_color(1.0, 0.0, 0.0);
        sphere
    })
    .collect();

    window.set_light(Light::StickToCamera);

    while window.render() {
        draw_box(
            &Point3::new(0., 1., 0.),
            max_x as f32,
            max_y as f32,
            max_z as f32,
            &mut window,
        );
        sphere_positions_x =
            rustplay::iter_integrate(time_step, &sphere_velocities_x[..], &sphere_positions_x[..])
                .collect();
        rustplay::iter_apply_bounds(
            0.,
            max_x,
            &mut sphere_velocities_x[..],
            &mut sphere_positions_x[..],
        );
        sphere_positions_y =
            rustplay::iter_integrate(time_step, &sphere_velocities_y[..], &sphere_positions_y[..])
                .collect();
        rustplay::iter_apply_bounds(
            0.,
            max_x,
            &mut sphere_velocities_y[..],
            &mut sphere_positions_y[..],
        );
        sphere_positions_z =
            rustplay::iter_integrate(time_step, &sphere_velocities_z[..], &sphere_positions_z[..])
                .collect();
        rustplay::iter_apply_bounds(
            0.,
            max_x,
            &mut sphere_velocities_z[..],
            &mut sphere_positions_z[..],
        );
        spheres = izip!(
            spheres,
            &sphere_positions_x,
            &sphere_positions_y,
            &sphere_positions_z,
        )
        .map(|(mut sphere, x, y, z)| {
            sphere.set_local_translation(Translation3::new(*x as f32, *y as f32, *z as f32));
            sphere
        })
        .collect();
        thread::sleep(draw_frequency);
    }
}
