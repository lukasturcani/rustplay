use rand::rngs::StdRng;
use rand::SeedableRng;
use rustplay::physics::spheres::{self, PhysicsConfig, XYZ};
use std::ptr;

fn main() {
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
    let max_velocity: f32 = 20.;
    let num_spheres = 350;
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

    unsafe {
        let flag = true;
        let flag = ptr::read_volatile(&flag);
        while flag {
            spheres::mix_take_time_step(&physics_config, &mut physics_spheres);
        }
    }
    println!("{:#?}", physics_spheres);
}
