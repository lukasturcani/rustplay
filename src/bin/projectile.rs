use kiss3d::nalgebra::{Point3, Translation3};

use kiss3d::light::Light;
use kiss3d::window::Window;

use clap::Parser;
use env_logger::Builder;
use itertools::izip;
use log::info;
use log::LevelFilter;
use rustplay::graphics;
use rustplay::physics::spheres::{self, PhysicsConfig, PhysicsData, XYZ};
use serde::Deserialize;
use serde_dhall;
use std::time::Instant;

#[derive(Parser, Debug)]
struct Args {
    /// Path to the config file.
    #[clap(value_parser)]
    config: String,
}

#[derive(Debug, Deserialize)]
enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

#[derive(Debug, Deserialize)]
struct LatticeConfig {
    dimensions: XYZ<usize>,
    offset: XYZ<f32>,
}

#[derive(Debug, Deserialize)]
struct Config {
    physics: PhysicsConfig,
    lattice: LatticeConfig,
    log: LogLevel,
}

fn log_level(level: &LogLevel) -> LevelFilter {
    match level {
        LogLevel::Error => LevelFilter::Error,
        LogLevel::Warn => LevelFilter::Warn,
        LogLevel::Info => LevelFilter::Info,
        LogLevel::Debug => LevelFilter::Debug,
        LogLevel::Trace => LevelFilter::Trace,
    }
}

fn add_projectile(physics_spheres: &mut PhysicsData) {
    physics_spheres.positions_x.push(5.);
    physics_spheres.positions_y.push(17.);
    physics_spheres.positions_z.push(17.);
    physics_spheres.velocities_x.push(10.);
    physics_spheres.velocities_y.push(0.);
    physics_spheres.velocities_z.push(0.);
}

fn main() -> Result<(), serde_dhall::Error> {
    let args = Args::parse();
    let config: Config = serde_dhall::from_file(args.config).parse()?;
    Builder::from_default_env()
        .filter(None, log_level(&config.log))
        .init();
    info!("Running with config {:#?}", config);
    let mut physics_spheres = spheres::get_lattice(
        (
            config.lattice.dimensions.x,
            config.lattice.dimensions.y,
            config.lattice.dimensions.z,
        ),
        (
            config.lattice.offset.x,
            config.lattice.offset.y,
            config.lattice.offset.z,
        ),
    );
    add_projectile(&mut physics_spheres);
    let mut window = Window::new("rustplay");
    let mut rendered_spheres: Vec<_> = izip!(
        &physics_spheres.positions_x,
        &physics_spheres.positions_y,
        &physics_spheres.positions_z,
    )
    .map(|(x, y, z)| {
        let mut sphere = window.add_sphere(config.physics.sphere_radius);
        sphere.set_local_translation(Translation3::new(*x, *y, *z));
        sphere.set_color(1.0, 0.0, 0.0);
        sphere
    })
    .collect();
    if let Some(sphere) = rendered_spheres.last_mut() {
        sphere.set_color(0., 1., 0.);
    }

    window.set_light(Light::StickToCamera);

    while window.render() {
        graphics::draw_box(
            &Point3::new(0., 1., 0.),
            config.physics.box_size.x,
            config.physics.box_size.y,
            config.physics.box_size.z,
            &mut window,
        );
        let physics_start = Instant::now();
        let mut simulated_time = 0.;
        while simulated_time < 0.016 {
            simulated_time +=
                spheres::exact_time_step(&config.physics, simulated_time, &mut physics_spheres);
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
    Ok(())
}
