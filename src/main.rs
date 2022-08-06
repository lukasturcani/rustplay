use kiss3d::nalgebra::{Point3, Translation3};

use kiss3d::light::Light;
use kiss3d::window::Window;

use itertools::izip;
use log::info;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rustplay::graphics;
use rustplay::physics::spheres::{self, PhysicsConfig, XYZ};
use std::time::Instant;

fn main() {
    env_logger::init();

    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        box_size: XYZ {
            x: 30.,
            y: 30.,
            z: 30.,
        },
        sphere_radius: 1.,
    };
    let max_velocity = 50.;
    let num_spheres = 50;
    let mut physics_spheres = spheres::get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.box_size.x,
        physics_config.box_size.y,
        physics_config.box_size.z,
        max_velocity,
        max_velocity,
        max_velocity,
    );
    let mut window = Window::new("rustplay");
    let mut rendered_spheres: Vec<_> = izip!(
        &physics_spheres.positions_x,
        &physics_spheres.positions_y,
        &physics_spheres.positions_z,
    )
    .map(|(x, y, z)| {
        let mut sphere = window.add_sphere(physics_config.sphere_radius);
        sphere.set_local_translation(Translation3::new(*x, *y, *z));
        sphere.set_color(1.0, 0.0, 0.0);
        sphere
    })
    .collect();

    window.set_light(Light::StickToCamera);

    while window.render() {
        graphics::draw_box(
            &Point3::new(0., 1., 0.),
            physics_config.box_size.x,
            physics_config.box_size.y,
            physics_config.box_size.z,
            &mut window,
        );
        let physics_start = Instant::now();
        let mut simulated_time = 0.;
        while simulated_time < 0.016 {
            simulated_time +=
                spheres::exact_time_step(&physics_config, simulated_time, &mut physics_spheres);
        }
        info!("Physics steps took {:#?}.", Instant::now() - physics_start);
        rendered_spheres = izip!(
            rendered_spheres,
            &physics_spheres.positions_x,
            &physics_spheres.positions_y,
            &physics_spheres.positions_z,
        )
        .map(|(mut sphere, x, y, z)| {
            sphere.set_local_translation(Translation3::new(*x, *y, *z));
            sphere
        })
        .collect();
    }
}
