module TestPath

using Modia
using Modia.Measurements
@usingModiaPlot

const ptp_path = PTP_path(["angle1", "angle2", "angle3"],
                          positions = [0.0 2.0 3.0;  # angle1=0.0, angle2=2.0, angle3=3.0
                                       0.5 3.0 4.0;
                                       0.8 1.5 0.3;
                                       0.2 1.5 0.8],
                          startTime = (0.1/60)u"minute",
                          v_max = 2*ones(3),
                          a_max = 3*ones(3))
angles = zeros(3)
getPosition!(ptp_path, 0.5, angles)   # angles = [0.12, 2.24, 3.24]
path = getPath(ptp_path)
showInfo(path)

plot(path, [("angle1", "angle2", "angle3"),
            ("der(angle1)", "der(angle2)", "der(angle3)"),
            ("der2(angle1)", "der2(angle2)", "der2(angle3)")], figure=1)

plotPath(ptp_path,plot,onlyPositions=false, heading="ptp_path", figure=2)


nom(v) = measurement(v)
err(v) = v ± (0.1*v)
const ptp_path2 = PTP_path{Measurement{Float64}}(
                          ["angle1", "angle2", "angle3"],
                           positions = [nom(0.0) nom(2.0) nom(3.0);  # angle1=0.0, angle2=2.0, angle3=3.0
                                        err(0.5) err(3.0) err(4.0);
                                        err(0.8) err(1.5) err(0.3);
                                        err(0.2) err(1.5) err(0.8)],
                           startTime = nom(0.1),
                           v_max = fill(measurement(2.0),3),
                           a_max = fill(measurement(3.0),3))
angles = zeros(Measurement{Float64}, 3)
getPosition!(ptp_path2, nom(0.5), angles)   # angles = [0.12, 2.24, 3.24]
path2 = getPath(ptp_path2)
showInfo(path2)

plot(path2, [("angle1", "angle2", "angle3"),
             ("der(angle1)", "der(angle2)", "der(angle3)"),
             ("der2(angle1)", "der2(angle2)", "der2(angle3)")], figure=3)
             

plotPath(ptp_path2,plot,onlyPositions=false, heading="ptp_path2", figure=4)
end