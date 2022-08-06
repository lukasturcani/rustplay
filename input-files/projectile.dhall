let XYZ = \(a : Type) -> { x : a, y : a, z : a }

let PhysicsConfig =
      { time_step : Double, box_size : XYZ Double, sphere_radius : Double }

let LatticeConfig = { dimensions : XYZ Natural, offset : XYZ Double }

let LogLevel = < Error | Warn | Info | Debug | Trace >

let config
    : { physics : PhysicsConfig, lattice : LatticeConfig, log : LogLevel }
    = { physics =
        { time_step = 0.0016
        , box_size = { x = 300.0, y = 300.0, z = 300.0 }
        , sphere_radius = 1.0
        }
      , lattice =
        { dimensions = { x = 2, y = 2, z = 2 }
        , offset = { x = 15.0, y = 15.0, z = 15.0 }
        }
      , log = LogLevel.Info
      }

in  config
