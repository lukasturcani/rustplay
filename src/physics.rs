use itertools::izip;

pub struct PhysicsConfig {
    pub time_step: f64,
    pub max_x: f64,
    pub max_y: f64,
    pub max_z: f64,
    pub sphere_radius: f64,
}

pub struct PhysicsData {
    pub positions_x: Vec<f64>,
    pub positions_y: Vec<f64>,
    pub positions_z: Vec<f64>,
    pub velocities_x: Vec<f64>,
    pub velocities_y: Vec<f64>,
    pub velocities_z: Vec<f64>,
}

pub fn get_random_physics_data(
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

struct SphereCollision {
    time: f64,
}

struct SphereCollisionId {
    sphere1: usize,
    sphere2: usize,
}

fn quadratic_root(a: f64, b: f64, c: f64) -> f64 {
    let term1 = -b;
    let term2 = (b * b - 4f64 * a * c).sqrt();
    let term3 = a + a;
    let root1 = (term1 + term2) / term3;
    let root2 = (term1 - term2) / term3;
    if root1 < 0f64 {
        root1
    } else if root2 < 0f64 {
        root1
    } else if root1 < root2 {
        root1
    } else {
        root2
    }
}

fn sphere_collision(
    r1: f64,
    px1: f64,
    py1: f64,
    pz1: f64,
    r2: f64,
    px2: f64,
    py2: f64,
    pz2: f64,
) -> bool {
    let x_diff = px2 - px1;
    let y_diff = py2 - py1;
    let z_diff = pz2 - pz1;
    let distance_squared = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
    let collision_distance = r1 + r2;
    distance_squared < collision_distance * collision_distance
}

fn sphere_collision_response(
    nx: f64,
    ny: f64,
    nz: f64,
    vx1: &mut f64,
    vy1: &mut f64,
    vz1: &mut f64,
    vx2: &mut f64,
    vy2: &mut f64,
    vz2: &mut f64,
) {
    let projection_magnitude1 = nx * *vx1 + ny * *vy1 + nz * *vz1;
    let projection_x1 = projection_magnitude1 * nx;
    let projection_y1 = projection_magnitude1 * ny;
    let projection_z1 = projection_magnitude1 * nz;

    let projection_magnitude2 = nx * *vx2 + ny * *vy2 + nz * *vz2;
    let projection_x2 = projection_magnitude2 * nx;
    let projection_y2 = projection_magnitude2 * ny;
    let projection_z2 = projection_magnitude2 * nz;
    *vx1 = *vx1 + projection_x2 - projection_x1;
    *vy1 = *vy1 + projection_y2 - projection_y1;
    *vz1 = *vz1 + projection_z2 - projection_z1;
    *vx2 = *vx2 + projection_x1 - projection_x2;
    *vy2 = *vy2 + projection_y1 - projection_y2;
    *vz2 = *vz2 + projection_z1 - projection_z2;
}

pub fn iter_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) -> f64 {
    let positions_x: Vec<_> = iter_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    )
    .collect();
    let positions_y: Vec<_> = iter_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    )
    .collect();
    let positions_z: Vec<_> = iter_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    )
    .collect();
    let first_collision = izip!(
        &positions_x,
        &positions_y,
        &positions_z,
        &data.velocities_x,
        &data.velocities_y,
        &data.velocities_z
    )
    .enumerate()
    .map(|(i1, (px1, py1, pz1, _, _, _))| {
        izip!(
            &positions_x,
            &positions_y,
            &positions_z,
            &data.velocities_x,
            &data.velocities_y,
            &data.velocities_z
        )
        .enumerate()
        .skip(i1 + 1)
        .filter(|(_, (px2, py2, pz2, _, _, _))| {
            sphere_collision(
                config.sphere_radius,
                *px1,
                *py1,
                *pz1,
                config.sphere_radius,
                **px2,
                **py2,
                **pz2,
            )
        })
        .map(|(i2, _)| {
            let px_diff = data.positions_x[i2] - data.positions_x[i1];
            let py_diff = data.positions_y[i2] - data.positions_y[i1];
            let pz_diff = data.positions_z[i2] - data.positions_z[i1];
            let vx_diff = data.velocities_x[i2] - data.velocities_x[i1];
            let vy_diff = data.velocities_y[i2] - data.velocities_y[i1];
            let vz_diff = data.velocities_z[i2] - data.velocities_z[i1];
            let distance_squared = px_diff * px_diff + py_diff * py_diff + pz_diff * pz_diff;
            let relative_speed_squared = vx_diff * vx_diff + vy_diff * vy_diff + vz_diff * vz_diff;
            let collision_distance = config.sphere_radius + config.sphere_radius;
            let half_b = vx_diff * px_diff + vy_diff * py_diff + vz_diff * pz_diff;
            let c = distance_squared - collision_distance * collision_distance;

            let time = quadratic_root(relative_speed_squared, half_b + half_b, c);
            if time.is_nan() {
                None
            } else {
                Some((
                    SphereCollisionId {
                        sphere1: i1,
                        sphere2: i2,
                    },
                    SphereCollision { time },
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

    let time = match first_collision {
        Some(collision) => {
            data.positions_x = iter_integrate(
                collision.1.time,
                &data.velocities_x[..],
                &data.positions_x[..],
            )
            .collect();
            data.positions_y = iter_integrate(
                collision.1.time,
                &data.velocities_y[..],
                &data.positions_y[..],
            )
            .collect();
            data.positions_z = iter_integrate(
                collision.1.time,
                &data.velocities_z[..],
                &data.positions_z[..],
            )
            .collect();
            let displacement = (
                data.positions_x[collision.0.sphere2] - data.positions_x[collision.0.sphere1],
                data.positions_y[collision.0.sphere2] - data.positions_y[collision.0.sphere1],
                data.positions_z[collision.0.sphere2] - data.positions_z[collision.0.sphere1],
            );
            let distance = (displacement.0 * displacement.0
                + displacement.1 * displacement.1
                + displacement.2 * displacement.2)
                .sqrt();
            let collision_normal = (
                displacement.0 / distance,
                displacement.1 / distance,
                displacement.2 / distance,
            );
            let (left_x, right_x) = data.velocities_x.split_at_mut(collision.0.sphere1 + 1);
            let (left_y, right_y) = data.velocities_y.split_at_mut(collision.0.sphere1 + 1);
            let (left_z, right_z) = data.velocities_z.split_at_mut(collision.0.sphere1 + 1);
            sphere_collision_response(
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
            collision.1.time
        }
        None => {
            data.positions_x = positions_x;
            data.positions_y = positions_y;
            data.positions_z = positions_z;
            config.time_step
        }
    };

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
    time
}

pub fn loop_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) {
    data.positions_x = loop_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    );
    loop_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );
    data.positions_y = loop_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    );
    loop_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    data.positions_z = loop_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    );
    loop_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
}

fn loop_integrate(time_step: f64, positions: &[f64], velocities: &[f64]) -> Vec<f64> {
    let mut output = Vec::with_capacity(positions.len());
    for i in 0..positions.len() {
        output.push(positions[i] + time_step * velocities[i]);
    }
    output
}

fn iter_integrate<'a>(
    time_step: f64,
    velocities: &'a [f64],
    positions: &'a [f64],
) -> impl Iterator<Item = f64> + 'a {
    let time_step_ = time_step.clone();
    positions
        .iter()
        .zip(velocities)
        .map(move |(position, velocity)| position + time_step_ * velocity)
}

fn loop_apply_bounds<'a>(
    min: f64,
    max: f64,
    velocities: &mut [f64],
    positions: &'a mut [f64],
) -> &'a mut [f64] {
    for i in 0..positions.len() {
        if positions[i] < min {
            positions[i] = min + min - positions[i];
            velocities[i] = -velocities[i];
        } else if positions[i] > max {
            positions[i] = max + max - positions[i];
            velocities[i] = -velocities[i];
        }
    }
    positions
}

fn iter_apply_bounds<'a>(
    min: f64,
    max: f64,
    velocities: &mut [f64],
    positions: &'a mut [f64],
) -> &'a mut [f64] {
    positions
        .iter_mut()
        .zip(velocities)
        .for_each(|(position, velocity)| {
            if *position < min {
                *position = &min + &min - *position;
                *velocity = -*velocity;
            } else if *position > max {
                *position = max + max - *position;
                *velocity = -*velocity;
            }
        });
    positions
}
