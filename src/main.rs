use kiss3d::nalgebra::{Point3, Translation3};

use kiss3d::light::Light;
use kiss3d::window::Window;

use itertools::izip;
use rand::rngs::StdRng;
use rand::SeedableRng;
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

struct PhysicsData {
    positions_x: Vec<f64>,
    positions_y: Vec<f64>,
    positions_z: Vec<f64>,
    velocities_x: Vec<f64>,
    velocities_y: Vec<f64>,
    velocities_z: Vec<f64>,
}

fn get_random_physics_data(
    generator: &mut impl rand::Rng,
    num_objects: u64,
    max_position_x: f64,
    max_position_y: f64,
    max_position_z: f64,
    max_velocity_x: f64,
    max_velocity_y: f64,
    max_velocity_z: f64,
) -> PhysicsData {
    PhysicsData {
        positions_x: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_x)
            .collect(),
        positions_y: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_y)
            .collect(),
        positions_z: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_z)
            .collect(),
        velocities_x: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_x)
            .collect(),
        velocities_y: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_y)
            .collect(),
        velocities_z: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_z)
            .collect(),
    }
}

struct PhysicsConfig {
    time_step: f64,
    max_x: f64,
    max_y: f64,
    max_z: f64,
}

fn take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) {
    data.positions_x = rustplay::iter_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    )
    .collect();
    rustplay::iter_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );
    data.positions_y = rustplay::iter_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    )
    .collect();
    rustplay::iter_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    data.positions_z = rustplay::iter_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    )
    .collect();
    rustplay::iter_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
}

fn main() {
    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 1.,
        max_x: 60.,
        max_y: 60.,
        max_z: 60.,
    };
    let max_velocity: f64 = 2.;
    let draw_frequency = time::Duration::from_millis(15);
    let num_spheres = 50;
    let mut physics_spheres = get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.max_x,
        physics_config.max_y,
        physics_config.max_z,
        max_velocity,
        max_velocity,
        max_velocity,
    );
    let mut window = Window::new("Kiss3d: cube");
    let mut rendered_spheres: Vec<_> = izip!(
        &physics_spheres.positions_x,
        &physics_spheres.positions_y,
        &physics_spheres.positions_z,
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
            physics_config.max_x as f32,
            physics_config.max_y as f32,
            physics_config.max_z as f32,
            &mut window,
        );
        take_time_step(&physics_config, &mut physics_spheres);
        rendered_spheres = izip!(
            rendered_spheres,
            &physics_spheres.positions_x,
            &physics_spheres.positions_y,
            &physics_spheres.positions_z,
        )
        .map(|(mut sphere, x, y, z)| {
            sphere.set_local_translation(Translation3::new(*x as f32, *y as f32, *z as f32));
            sphere
        })
        .collect();
        thread::sleep(draw_frequency);
    }
}
