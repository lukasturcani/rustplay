let XYZ = { x : Double, y : Double, z : Double }

let PhysicsConfig =
      { time_step : Double, box_size : XYZ, sphere_radius : Double }

let LatticeConfig = { dimensions : XYZ, offset : XYZ }

let LogLevel = < Error | Warn | Info | Debug | Trace >

let config
    : { physics : PhysicsConfig, lattice : LatticeConfig, log : LogLevel }
    = { physics =
        { time_step = 0.0016
        , box_size = { x = 300.0, y = 300.0, z = 300.0 }
        , sphere_radius = 1.0
        }
      , lattice =
        { dimensions = { x = 2.0, y = 2.0, z = 2.0 }
        , offset = { x = 15.0, y = 15.0, z = 15.0 }
        }
      , log = LogLevel.Info
      }

in  config
