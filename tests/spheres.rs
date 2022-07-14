use itertools::izip;
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Translation3};
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use rustplay::graphics;
use rustplay::physics::{self, PhysicsConfig, PhysicsData};

fn get_window() -> Window {
    let mut window = Window::new("rustplay");
    window.set_light(Light::StickToCamera);
    window
}

fn add_spheres<'a>(
    physics_config: &'a PhysicsConfig,
    physics_spheres: &'a PhysicsData,
    window: &'a mut Window,
) -> impl Iterator<Item = SceneNode> + 'a {
    izip!(
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
}

fn physics_loop(
    simulation_time: f32,
    window: &mut Window,
    physics_config: &PhysicsConfig,
    physics_spheres: &mut PhysicsData,
    mut rendered_spheres: Vec<SceneNode>,
) {
    let mut total_simulated_time = 0.;
    while window.render() && total_simulated_time < simulation_time {
        graphics::draw_box(
            &Point3::new(0., 1., 0.),
            physics_config.max_x,
            physics_config.max_y,
            physics_config.max_z,
            window,
        );
        let mut simulated_time = 0.;
        while simulated_time < 0.016 {
            simulated_time +=
                physics::iter_take_time_step(physics_config, simulated_time, physics_spheres);
        }
        total_simulated_time += simulated_time;
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

#[test]
fn collision() {
    let physics_config = PhysicsConfig {
        time_step: 0.0016,
        max_x: 30.,
        max_y: 30.,
        max_z: 30.,
        sphere_radius: 1.,
    };
    let mut physics_spheres = PhysicsData {
        positions_x: vec![50., 100.],
        positions_y: vec![0., 0.],
        positions_z: vec![0., 0.],
        velocities_x: vec![10., -10.],
        velocities_y: vec![0., 0.],
        velocities_z: vec![0., 0.],
    };
    let mut window = get_window();
    let rendered_spheres: Vec<_> =
        add_spheres(&physics_config, &physics_spheres, &mut window).collect();
    physics_loop(
        25.,
        &mut window,
        &physics_config,
        &mut physics_spheres,
        rendered_spheres,
    );
}
