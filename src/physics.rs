use core::arch::x86_64;
use itertools::izip;

pub struct PhysicsConfig {
    pub time_step: f32,
    pub max_x: f32,
    pub max_y: f32,
    pub max_z: f32,
    pub sphere_radius: f32,
}

#[derive(Debug)]
pub struct PhysicsData {
    pub positions_x: Vec<f32>,
    pub positions_y: Vec<f32>,
    pub positions_z: Vec<f32>,
    pub velocities_x: Vec<f32>,
    pub velocities_y: Vec<f32>,
    pub velocities_z: Vec<f32>,
}

pub fn get_random_physics_data(
    generator: &mut impl rand::Rng,
    num_objects: u64,
    max_position_x: f32,
    max_position_y: f32,
    max_position_z: f32,
    max_velocity_x: f32,
    max_velocity_y: f32,
    max_velocity_z: f32,
) -> PhysicsData {
    PhysicsData {
        positions_x: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_position_x)
            .collect(),
        positions_y: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_position_y)
            .collect(),
        positions_z: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_position_z)
            .collect(),
        velocities_x: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_velocity_x)
            .collect(),
        velocities_y: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_velocity_y)
            .collect(),
        velocities_z: (0..num_objects)
            .map(|_| generator.gen::<f32>() * max_velocity_z)
            .collect(),
    }
}

struct SphereCollision {
    time: f32,
}

struct SphereCollisionId {
    sphere1: usize,
    sphere2: usize,
}

fn quadratic_root(a: f32, b: f32, c: f32) -> f32 {
    let term1 = -b;
    let term2 = (b * b - 4f32 * a * c).sqrt();
    let term3 = a + a;
    let root1 = (term1 + term2) / term3;
    let root2 = (term1 - term2) / term3;
    if root1 < 0f32 {
        root2
    } else if root2 < 0f32 {
        root1
    } else if root1 < root2 {
        root1
    } else {
        root2
    }
}

fn vector_sub(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 - x2, y1 - y2, z1 - z2)
}

fn vector_sub_((x1, y1, z1): (f32, f32, f32), (x2, y2, z2): (f32, f32, f32)) -> (f32, f32, f32) {
    vector_sub(x1, y1, z1, x2, y2, z2)
}

fn vector_add(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 + x2, y1 + y2, z1 + z2)
}

fn vector_add_((x1, y1, z1): (f32, f32, f32), (x2, y2, z2): (f32, f32, f32)) -> (f32, f32, f32) {
    vector_add(x1, y1, z1, x2, y2, z2)
}

fn vector_mul(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 * x2, y1 * y2, z1 * z2)
}

fn vector_scalar_mul(scalar: f32, x: f32, y: f32, z: f32) -> (f32, f32, f32) {
    (scalar * x, scalar * y, scalar * z)
}

fn vector_sum_((x, y, z): (f32, f32, f32)) -> f32 {
    x + y + z
}

fn vector_dot(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> f32 {
    vector_sum_(vector_mul(x1, y1, z1, x2, y2, z2))
}

fn vector_dot_((x1, y1, z1): (f32, f32, f32), (x2, y2, z2): (f32, f32, f32)) -> f32 {
    vector_dot(x1, y1, z1, x2, y2, z2)
}

fn is_sphere_collision(
    r1: f32,
    px1: f32,
    py1: f32,
    pz1: f32,
    r2: f32,
    px2: f32,
    py2: f32,
    pz2: f32,
) -> bool {
    let diff = vector_sub(px2, py2, pz2, px1, py1, pz1);
    let distance_squared = vector_dot_(diff, diff);
    let collision_distance = r1 + r2;
    distance_squared < collision_distance * collision_distance
}

fn sphere_collision_time(
    r1: f32,
    px1: f32,
    py1: f32,
    pz1: f32,
    vx1: f32,
    vy1: f32,
    vz1: f32,
    r2: f32,
    px2: f32,
    py2: f32,
    pz2: f32,
    vx2: f32,
    vy2: f32,
    vz2: f32,
) -> f32 {
    let p_diff = vector_sub(px2, py2, pz2, px1, py1, pz1);
    let v_diff = vector_sub(vx2, vy2, vz2, vx1, vy1, vz1);
    let distance_squared = vector_dot_(p_diff, p_diff);
    let relative_speed_squared = vector_dot_(v_diff, v_diff);
    let collision_distance = r1 + r2;
    let half_b = vector_dot_(p_diff, v_diff);
    let c = distance_squared - collision_distance * collision_distance;
    quadratic_root(relative_speed_squared, half_b + half_b, c)
}

fn sphere_collision_response(
    nx: f32,
    ny: f32,
    nz: f32,
    vx1: &mut f32,
    vy1: &mut f32,
    vz1: &mut f32,
    vx2: &mut f32,
    vy2: &mut f32,
    vz2: &mut f32,
) {
    let projection_magnitude1 = vector_dot(nx, ny, nz, *vx1, *vy1, *vz1);
    let projection1 = vector_scalar_mul(projection_magnitude1, nx, ny, nz);

    let projection_magnitude2 = vector_dot(nx, ny, nz, *vx2, *vy2, *vz2);
    let projection2 = vector_scalar_mul(projection_magnitude2, nx, ny, nz);
    (*vx1, *vy1, *vz1) = vector_add_((*vx1, *vy1, *vz1), vector_sub_(projection2, projection1));
    (*vx2, *vy2, *vz2) = vector_add_((*vx2, *vy2, *vz2), vector_sub_(projection1, projection2));
}

pub fn mix_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) -> f32 {
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

    let mut first_collision = (f32::INFINITY, 0, 0);
    for i1 in 0..positions_x.len() - 1 {
        for i2 in i1 + 1..positions_x.len() {
            if is_sphere_collision(
                config.sphere_radius,
                positions_x[i1],
                positions_y[i1],
                positions_z[i1],
                config.sphere_radius,
                positions_x[i2],
                positions_y[i2],
                positions_z[i2],
            ) {
                let time = sphere_collision_time(
                    config.sphere_radius,
                    data.positions_x[i1],
                    data.positions_y[i1],
                    data.positions_z[i1],
                    data.velocities_x[i1],
                    data.velocities_y[i1],
                    data.velocities_z[i1],
                    config.sphere_radius,
                    data.positions_x[i2],
                    data.positions_y[i2],
                    data.positions_z[i2],
                    data.velocities_x[i2],
                    data.velocities_y[i2],
                    data.velocities_z[i2],
                );
                if time.is_normal() && time < first_collision.0 {
                    first_collision = (time, i1, i2);
                }
            }
        }
    }
    let simulated_time = if first_collision.0.is_normal() {
        data.positions_x = iter_integrate(
            first_collision.0,
            &data.velocities_x[..],
            &data.positions_x[..],
        )
        .collect();
        data.positions_y = iter_integrate(
            first_collision.0,
            &data.velocities_y[..],
            &data.positions_y[..],
        )
        .collect();
        data.positions_z = iter_integrate(
            first_collision.0,
            &data.velocities_z[..],
            &data.positions_z[..],
        )
        .collect();
        let displacement = (
            data.positions_x[first_collision.2] - data.positions_x[first_collision.1],
            data.positions_y[first_collision.2] - data.positions_y[first_collision.1],
            data.positions_z[first_collision.2] - data.positions_z[first_collision.1],
        );
        let distance = vector_dot_(displacement, displacement).sqrt();
        let collision_normal = (
            displacement.0 / distance,
            displacement.1 / distance,
            displacement.2 / distance,
        );
        let (left_x, right_x) = data.velocities_x.split_at_mut(first_collision.1 + 1);
        let (left_y, right_y) = data.velocities_y.split_at_mut(first_collision.1 + 1);
        let (left_z, right_z) = data.velocities_z.split_at_mut(first_collision.1 + 1);
        sphere_collision_response(
            collision_normal.0,
            collision_normal.1,
            collision_normal.2,
            &mut right_x[first_collision.2 - left_x.len()],
            &mut right_y[first_collision.2 - left_y.len()],
            &mut right_z[first_collision.2 - left_z.len()],
            &mut left_x[first_collision.1],
            &mut left_y[first_collision.1],
            &mut left_z[first_collision.1],
        );
        first_collision.0
    } else {
        config.time_step
    };

    loop_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );

    loop_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    loop_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
    simulated_time
}

pub fn iter_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) -> f32 {
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
        &data.positions_x,
        &data.positions_y,
        &data.positions_z,
        &data.velocities_x,
        &data.velocities_y,
        &data.velocities_z,
        &positions_x,
        &positions_y,
        &positions_z,
    )
    .enumerate()
    .map(
        |(i1, (p0x1, p0y1, p0z1, vx1, vy1, vz1, p1x1, p1y1, p1z1))| {
            izip!(
                &data.positions_x,
                &data.positions_y,
                &data.positions_z,
                &data.velocities_x,
                &data.velocities_y,
                &data.velocities_z,
                &positions_x,
                &positions_y,
                &positions_z,
            )
            .enumerate()
            .skip(i1 + 1)
            .filter(|(_, (_, _, _, _, _, _, p1x2, p1y2, p1z2))| {
                is_sphere_collision(
                    config.sphere_radius,
                    *p1x1,
                    *p1y1,
                    *p1z1,
                    config.sphere_radius,
                    **p1x2,
                    **p1y2,
                    **p1z2,
                )
            })
            .map(|(i2, (p0x2, p0y2, p0z2, vx2, vy2, vz2, _, _, _))| {
                let time = sphere_collision_time(
                    config.sphere_radius,
                    *p0x1,
                    *p0y1,
                    *p0z1,
                    *vx1,
                    *vy1,
                    *vz1,
                    config.sphere_radius,
                    *p0x2,
                    *p0y2,
                    *p0z2,
                    *vx2,
                    *vy2,
                    *vz2,
                );
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
        },
    )
    .filter(Option::is_some)
    .map(Option::unwrap)
    .reduce(|first_collision, second_collision| {
        if first_collision.1.time < second_collision.1.time {
            first_collision
        } else {
            second_collision
        }
    });

    let simulated_time = match first_collision {
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
            let distance = vector_dot_(displacement, displacement).sqrt();
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
    simulated_time
}

pub fn loop_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) -> f32 {
    let positions_x = loop_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    );
    let positions_y = loop_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    );
    let positions_z = loop_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    );

    let mut first_collision = (f32::INFINITY, 0, 0);
    for i1 in 0..positions_x.len() - 1 {
        for i2 in i1 + 1..positions_x.len() {
            if is_sphere_collision(
                config.sphere_radius,
                positions_x[i1],
                positions_y[i1],
                positions_z[i1],
                config.sphere_radius,
                positions_x[i2],
                positions_y[i2],
                positions_z[i2],
            ) {
                let time = sphere_collision_time(
                    config.sphere_radius,
                    data.positions_x[i1],
                    data.positions_y[i1],
                    data.positions_z[i1],
                    data.velocities_x[i1],
                    data.velocities_y[i1],
                    data.velocities_z[i1],
                    config.sphere_radius,
                    data.positions_x[i2],
                    data.positions_y[i2],
                    data.positions_z[i2],
                    data.velocities_x[i2],
                    data.velocities_y[i2],
                    data.velocities_z[i2],
                );
                if time.is_normal() && time < first_collision.0 {
                    first_collision = (time, i1, i2);
                }
            }
        }
    }
    let simulated_time = if first_collision.0.is_normal() {
        data.positions_x = loop_integrate(
            first_collision.0,
            &data.velocities_x[..],
            &data.positions_x[..],
        );
        data.positions_y = loop_integrate(
            first_collision.0,
            &data.velocities_y[..],
            &data.positions_y[..],
        );
        data.positions_z = loop_integrate(
            first_collision.0,
            &data.velocities_z[..],
            &data.positions_z[..],
        );
        let displacement = (
            data.positions_x[first_collision.2] - data.positions_x[first_collision.1],
            data.positions_y[first_collision.2] - data.positions_y[first_collision.1],
            data.positions_z[first_collision.2] - data.positions_z[first_collision.1],
        );
        let distance = vector_dot_(displacement, displacement).sqrt();
        let collision_normal = (
            displacement.0 / distance,
            displacement.1 / distance,
            displacement.2 / distance,
        );
        let (left_x, right_x) = data.velocities_x.split_at_mut(first_collision.1 + 1);
        let (left_y, right_y) = data.velocities_y.split_at_mut(first_collision.1 + 1);
        let (left_z, right_z) = data.velocities_z.split_at_mut(first_collision.1 + 1);
        sphere_collision_response(
            collision_normal.0,
            collision_normal.1,
            collision_normal.2,
            &mut right_x[first_collision.2 - left_x.len()],
            &mut right_y[first_collision.2 - left_y.len()],
            &mut right_z[first_collision.2 - left_z.len()],
            &mut left_x[first_collision.1],
            &mut left_y[first_collision.1],
            &mut left_z[first_collision.1],
        );
        first_collision.0
    } else {
        config.time_step
    };

    loop_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );

    loop_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    loop_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
    simulated_time
}

fn loop_integrate(time_step: f32, velocities: &[f32], positions: &[f32]) -> Vec<f32> {
    let mut output = Vec::with_capacity(positions.len());
    for i in 0..positions.len() {
        output.push(positions[i] + time_step * velocities[i]);
    }
    output
}

fn iter_integrate<'a>(
    time_step: f32,
    velocities: &'a [f32],
    positions: &'a [f32],
) -> impl Iterator<Item = f32> + 'a {
    positions
        .iter()
        .zip(velocities)
        .map(move |(position, velocity)| position + time_step * velocity)
}

fn _simd_integrate(time_step: f32, velocities: &[f32], positions: &[f32]) -> Vec<f32> {
    let num_items = positions.len();
    let mut output = vec![0f32; num_items];
    unsafe {
        let time_step_ = x86_64::_mm256_set_ps(
            time_step, time_step, time_step, time_step, time_step, time_step, time_step, time_step,
        );
        for i in (0..num_items).step_by(8) {
            let position = x86_64::_mm256_loadu_ps(&positions[i]);
            let velocity = x86_64::_mm256_loadu_ps(&velocities[i]);
            let displacement = x86_64::_mm256_mul_ps(velocity, time_step_);
            let new_position = x86_64::_mm256_add_ps(position, displacement);
            x86_64::_mm256_storeu_ps(&mut output[i], new_position);
        }
        let left_over = num_items % 8;
        for i in (num_items - left_over)..num_items {
            output[i] = positions[i] + time_step * velocities[i];
        }
    }
    output
}

fn loop_apply_bounds<'a>(
    min: f32,
    max: f32,
    velocities: &mut [f32],
    positions: &'a mut [f32],
) -> &'a mut [f32] {
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
    min: f32,
    max: f32,
    velocities: &mut [f32],
    positions: &'a mut [f32],
) -> &'a mut [f32] {
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

#[cfg(test)]
mod tests {
    #[test]
    fn vector_mul() {
        assert_eq!(super::vector_mul(1., 2., 3., 0., 0., 0.), (0., 0., 0.));
        assert_eq!(super::vector_mul(1., 2., 3., 1., 1., 1.), (1., 2., 3.));
    }
    #[test]
    fn vector_dot() {
        assert_eq!(super::vector_dot(1., 2., 3., 1., 1., 1.), 6.);
        assert_eq!(super::vector_dot(1., 2., 3., 0., 0., 0.), 0.);
        assert_eq!(super::vector_dot(1., 2., 3., 1., 2., 3.), 14.);
        assert_eq!(super::vector_dot(1., 2., 3., 4., 5., 6.), 32.);
    }
    #[test]
    fn vector_dot_() {
        assert_eq!(super::vector_dot_((1., 2., 3.), (1., 1., 1.)), 6.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (0., 0., 0.)), 0.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (1., 2., 3.)), 14.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (4., 5., 6.)), 32.);
    }
    #[test]
    fn vector_sub() {
        assert_eq!(super::vector_sub(1., 2., 3., 0., 0., 0.), (1., 2., 3.));
        assert_eq!(super::vector_sub(1., 2., 3., 1., 2., 3.), (0., 0., 0.));
    }
    #[test]
    fn vector_sub_() {
        assert_eq!(super::vector_sub_((1., 2., 3.), (0., 0., 0.)), (1., 2., 3.));
        assert_eq!(super::vector_sub_((1., 2., 3.), (1., 2., 3.)), (0., 0., 0.));
    }
    #[test]
    fn vector_add() {
        assert_eq!(super::vector_add(1., 2., 3., 0., 0., 0.), (1., 2., 3.));
        assert_eq!(super::vector_add(1., 2., 3., 4., 5., 6.), (5., 7., 9.));
    }
    #[test]
    fn vector_add_() {
        assert_eq!(super::vector_add_((1., 2., 3.), (0., 0., 0.)), (1., 2., 3.));
        assert_eq!(super::vector_add_((1., 2., 3.), (4., 5., 6.)), (5., 7., 9.));
    }
    #[test]
    fn vector_sum_() {
        assert_eq!(super::vector_sum_((1., 2., 3.)), 6.);
    }
    #[test]
    fn vector_scalar_mul() {
        assert_eq!(super::vector_scalar_mul(2., 1., 2., 3.), (2., 4., 6.));
    }
    #[test]
    fn iter_apply_bounds() {
        let min = 0.;
        let max = 50.;
        let mut velocities = vec![1., 2., 3., 4.];
        let mut positions = vec![-1., 40., 53., 21.];
        super::iter_apply_bounds(min, max, &mut velocities, &mut positions);
        assert_eq!(velocities, vec![-1., 2., -3., 4.]);
        assert_eq!(positions, vec![1., 40., 47., 21.]);
    }
    #[test]
    fn sphere_collision_response() {
        let (nx, ny, nz) = (1., 0., 0.);
        let (mut vx1, mut vy1, mut vz1) = (1., 0., 0.);
        let (mut vx2, mut vy2, mut vz2) = (-1., 0., 0.);
        super::sphere_collision_response(
            nx, ny, nz, &mut vx1, &mut vy1, &mut vz1, &mut vx2, &mut vy2, &mut vz2,
        );
        assert_eq!((vx1, vy1, vz1), (-1., 0., 0.));
    }
}

#[cfg(feature = "bench")]
pub mod bench {
    pub fn iter_integrate<'a>(
        time_step: f32,
        velocities: &'a [f32],
        positions: &'a [f32],
    ) -> impl Iterator<Item = f32> + 'a {
        super::iter_integrate(time_step, &velocities, &positions)
    }

    pub fn iter_apply_bounds<'a>(
        min: f32,
        max: f32,
        velocities: &mut [f32],
        positions: &'a mut [f32],
    ) -> &'a mut [f32] {
        super::iter_apply_bounds(min, max, velocities, positions)
    }

    pub fn loop_integrate(time_step: f32, velocities: &[f32], positions: &[f32]) -> Vec<f32> {
        super::loop_integrate(time_step, &velocities, &positions)
    }

    pub fn loop_apply_bounds<'a>(
        min: f32,
        max: f32,
        velocities: &mut [f32],
        positions: &'a mut [f32],
    ) -> &'a mut [f32] {
        super::loop_apply_bounds(min, max, velocities, positions)
    }

    pub fn simd_integrate(time_step: f32, velocities: &[f32], positions: &[f32]) -> Vec<f32> {
        super::_simd_integrate(time_step, &velocities, &positions)
    }
}
