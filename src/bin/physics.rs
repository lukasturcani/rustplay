use rand::rngs::StdRng;
use rand::SeedableRng;
use rustplay::physics::{self, PhysicsConfig};
use std::ptr;

fn main() {
    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let max_velocity: f64 = 20.;
    let num_spheres = 350;
    let mut physics_spheres = physics::get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.max_x,
        physics_config.max_y,
        physics_config.max_z,
        max_velocity,
        max_velocity,
        max_velocity,
    );

    unsafe {
        let flag = true;
        let flag = ptr::read_volatile(&flag);
        while flag {
            physics::mix_take_time_step(&physics_config, &mut physics_spheres);
        }
    }
    println!("{:#?}", physics_spheres);
}
