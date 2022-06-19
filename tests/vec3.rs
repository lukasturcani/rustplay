use rustplay;

#[test]
fn add() {
    assert_eq!(
        rustplay::vec3::add(
            &rustplay::vec3::Vec3 {
                x: 1.,
                y: 2.,
                z: 3.
            },
            &rustplay::vec3::Vec3 {
                x: 4.,
                y: 5.,
                z: 6.
            },
        ),
        rustplay::vec3::Vec3 {
            x: 1. + 4.,
            y: 2. + 5.,
            z: 3. + 6.
        },
    );
}
