module TestPath

using TinyModia, ModiaPlot

const ptp_path = PTP_path(["angle1", "angle2", "angle3"],
                          positions = [0.0 2.0 3.0;  # angle1=0.0, angle2=2.0, angle3=3.0
                                       0.5 3.0 4.0;
                                       0.8 1.5 0.3;
                                       0.2 1.5 0.8],
                          startTime = 0.1,
                          v_max = 2*ones(3),
                          a_max = 3*ones(3))
angles = zeros(3)
getPosition!(ptp_path, 0.5, angles)   # angles = [0.12, 2.24, 3.24]
path = getPath(ptp_path)
printResultInfo(path)

plot(path, [("angle1", "angle2", "angle3"),
            ("der(angle1)", "der(angle2)", "der(angle3)"),
            ("der2(angle1)", "der2(angle2)", "der2(angle3)")])
            
end