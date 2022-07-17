use kiss3d::nalgebra::{Point3, Translation3};

use kiss3d::light::Light;
use kiss3d::window::Window;

use itertools::izip;
use log::info;
use rustplay::graphics;
use rustplay::physics::spheres::{self, PhysicsConfig, PhysicsData};
use std::time::Instant;

fn add_projectile(physics_spheres: &mut PhysicsData) {
    physics_spheres.positions_x.push(5.);
    physics_spheres.positions_y.push(17.);
    physics_spheres.positions_z.push(17.);
    physics_spheres.velocities_x.push(10.);
    physics_spheres.velocities_y.push(0.);
    physics_spheres.velocities_z.push(0.);
}

fn main() {
    env_logger::init();

    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let mut physics_spheres = spheres::get_lattice((5, 5, 5), (15., 15., 15.));
    add_projectile(&mut physics_spheres);
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
            physics_config.max_x,
            physics_config.max_y,
            physics_config.max_z,
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
