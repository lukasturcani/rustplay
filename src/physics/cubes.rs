// use super::common;
// use itertools::izip;

// pub struct PhysicsConfig {
//     time_step: f32,
//     cube_size: f32,
// }

// struct CubeCollisionId {
//     cube1: usize,
//     cube2: usize,
// }

// struct CubeCollision {
//     time: f32,
// }

#[derive(Debug)]
pub struct Quaternion {}

#[derive(Debug)]
pub struct PhysicsData {
    pub positions_x: Vec<f32>,
    pub positions_y: Vec<f32>,
    pub positions_z: Vec<f32>,
    pub velocities_x: Vec<f32>,
    pub velocities_y: Vec<f32>,
    pub velocities_z: Vec<f32>,
    pub orientations: Vec<Quaternion>,
}

/*
pub fn exact_time_step(config: &PhysicsConfig, simulated_time: f32, data: &mut PhysicsData) -> f32 {
    let mut time_to_simulate = config.time_step - simulated_time;
    if time_to_simulate < 0.0001 {
        time_to_simulate = 0.0001;
    }
    let time_to_simulate = time_to_simulate;

    let first_collision = izip!(
        &data.positions_x,
        &data.positions_y,
        &data.positions_z,
        &data.velocities_x,
        &data.velocities_y,
        &data.velocities_z,
    )
    .enumerate()
    .map(|(i1, (px1, py1, pz1, vx1, vy1, vz1))| {
        izip!(
            &data.positions_x,
            &data.positions_y,
            &data.positions_z,
            &data.velocities_x,
            &data.velocities_y,
            &data.velocities_z,
        )
        .enumerate()
        .skip(i1 + 1)
        .map(|(i2, (px2, py2, pz2, vx2, vy2, vz2))| {
            let time = cube_collision_time(
                config.cube_size,
                *px1,
                *py1,
                *pz1,
                *vx1,
                *vy1,
                *vz1,
                config.cube_size,
                *px2,
                *py2,
                *pz2,
                *vx2,
                *vy2,
                *vz2,
            );
            if time.is_nan() || time <= 0. {
                None
            } else {
                Some((
                    CubeCollisionId {
                        cube1: i1,
                        cube2: i2,
                    },
                    CubeCollision { time },
                ))
            }
        })
        .filter(Option::is_some)
        .map(Option::unwrap)
        .reduce(|first_collision, second_collision| {
            if first_collision.1.time < second_collision.1.time {
                first_collision
            } else {
                second_collision
            }
        })
    })
    .filter(Option::is_some)
    .map(Option::unwrap)
    .reduce(|first_collision, second_collision| {
        if first_collision.1.time < second_collision.1.time {
            first_collision
        } else {
            second_collision
        }
    });

    let simulated_time = match &first_collision {
        Some(collision) => {
            if collision.1.time < time_to_simulate {
                collision.1.time
            } else {
                time_to_simulate
            }
        }
        None => time_to_simulate,
    };
    data.positions_x = common::iter_integrate(
        simulated_time,
        &data.velocities_x[..],
        &data.positions_x[..],
    )
    .collect();
    data.positions_y = common::iter_integrate(
        simulated_time,
        &data.velocities_y[..],
        &data.positions_y[..],
    )
    .collect();
    data.positions_z = common::iter_integrate(
        simulated_time,
        &data.velocities_z[..],
        &data.positions_z[..],
    )
    .collect();

    if let Some(collision) = first_collision {
        let displacement = (
            data.positions_x[collision.0.sphere2] - data.positions_x[collision.0.sphere1],
            data.positions_y[collision.0.sphere2] - data.positions_y[collision.0.sphere1],
            data.positions_z[collision.0.sphere2] - data.positions_z[collision.0.sphere1],
        );
        let distance = common::vector_dot_(displacement, displacement).sqrt();
        let collision_normal = (
            displacement.0 / distance,
            displacement.1 / distance,
            displacement.2 / distance,
        );
        let (left_x, right_x) = data.velocities_x.split_at_mut(collision.0.sphere1 + 1);
        let (left_y, right_y) = data.velocities_y.split_at_mut(collision.0.sphere1 + 1);
        let (left_z, right_z) = data.velocities_z.split_at_mut(collision.0.sphere1 + 1);
        cube_collision_response(
            collision_normal.0,
            collision_normal.1,
            collision_normal.2,
            &mut right_x[collision.0.sphere2 - left_x.len()],
            &mut right_y[collision.0.sphere2 - left_y.len()],
            &mut right_z[collision.0.sphere2 - left_z.len()],
            &mut left_x[collision.0.sphere1],
            &mut left_y[collision.0.sphere1],
            &mut left_z[collision.0.sphere1],
        );
    }

    iter_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );
    iter_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    iter_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
    simulated_time
}
*/
