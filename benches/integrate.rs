use criterion::{criterion_group, criterion_main, Criterion};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rustplay::physics::{self, PhysicsConfig};

pub fn benchmark_iter_integrate(c: &mut Criterion) {
    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let max_velocity: f32 = 20.;
    let num_spheres = 350000;
    let physics_spheres = physics::get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.max_x,
        physics_config.max_y,
        physics_config.max_z,
        max_velocity,
        max_velocity,
        max_velocity,
    );
    c.bench_function("iter_integrate", |b| {
        b.iter(|| {
            let _result: Vec<_> = physics::bench::iter_integrate(
                physics_config.time_step,
                &physics_spheres.velocities_x,
                &physics_spheres.positions_x,
            )
            .collect();
        })
    });
}

pub fn benchmark_loop_integrate(c: &mut Criterion) {
    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let max_velocity: f32 = 20.;
    let num_spheres = 350000;
    let physics_spheres = physics::get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.max_x,
        physics_config.max_y,
        physics_config.max_z,
        max_velocity,
        max_velocity,
        max_velocity,
    );
    c.bench_function("loop_integrate", |b| {
        b.iter(|| {
            let _result = physics::bench::loop_integrate(
                physics_config.time_step,
                &physics_spheres.velocities_x,
                &physics_spheres.positions_x,
            );
        })
    });
}

pub fn benchmark_simd_integrate(c: &mut Criterion) {
    let mut generator = StdRng::seed_from_u64(12);
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let max_velocity: f32 = 20.;
    let num_spheres = 350000;
    let physics_spheres = physics::get_random_physics_data(
        &mut generator,
        num_spheres,
        physics_config.max_x,
        physics_config.max_y,
        physics_config.max_z,
        max_velocity,
        max_velocity,
        max_velocity,
    );
    c.bench_function("simd_integrate", |b| {
        b.iter(|| {
            let _result = physics::bench::simd_integrate(
                physics_config.time_step,
                &physics_spheres.velocities_x,
                &physics_spheres.positions_x,
            );
        })
    });
}

criterion_group!(
    benches,
    benchmark_iter_integrate,
    benchmark_loop_integrate,
    benchmark_simd_integrate,
);
criterion_main!(benches);
