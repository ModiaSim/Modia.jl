within ;
package ServoSystem
  "Servo system consisting of current and speed controlled DC motor + gear with elasticity and damping + load"

  model Gear "Gear with elasticity and damping"
    parameter Real ratio=105 "Getriebe-Uebersetzung";
    Modelica.Mechanics.Rotational.Components.IdealGear gear(
                                                 ratio=ratio, useSupport=true)
      annotation (Placement(transformation(extent={{-60,-10},{-40,10}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b
      annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=
           0)));
    Modelica.Mechanics.Rotational.Components.Spring
                    spring(c=5.84e5)
                                    annotation (Placement(transformation(extent=
           {{20,10},{40,30}}, rotation=0)));
    Modelica.Mechanics.Rotational.Components.Damper damper1(
                                                 d=500)
      annotation (Placement(transformation(extent={{20,-30},{40,-10}}, rotation=
           0)));
    Modelica.Mechanics.Rotational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-60,-49},{-40,-29}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Components.Damper damper2(
                                                 d=100)
      annotation (Placement(transformation(
        origin={-30,-20},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  equation
    connect(gear.flange_a, flange_a)
      annotation (Line(points={{-60,0},{-100,0}}, color={0,0,0}));
    connect(gear.flange_b, spring.flange_a)
      annotation (Line(points={{-40,0},{0,0},{0,20},{20,20}}, color={0,0,0}));
    connect(gear.flange_b, damper1.flange_a)
      annotation (Line(points={{-40,0},{0,0},{0,-20},{20,-20}}, color={0,0,0}));
    connect(spring.flange_b, flange_b) annotation (Line(points={{40,20},{60,20},
            {60,0},{100,0}}, color={0,0,0}));
    connect(damper1.flange_b, flange_b)
      annotation (Line(points={{40,-20},{60,-20},{60,0},{100,0}}, color={0,0,0}));
    connect(fixed.flange,   damper2.flange_b)
      annotation (Line(points={{-50,-39},{-50,-30},{-30,-30}}, color={0,0,0}));
    connect(gear.support,fixed.flange)
      annotation (Line(points={{-50,-10},{-50,-39}}, color={0,0,0}));
    connect(damper2.flange_a, gear.flange_b) annotation (Line(points={{-30,-10},{
            -30,0},{-40,0}},color={0,0,0}));
    annotation (
      Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics={
          Rectangle(
            extent={{-100,10},{-60,-10}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{60,10},{100,-10}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-40,60},{40,-60}},
            lineColor={0,0,0},
            pattern=LinePattern.Solid,
            lineThickness=0.25,
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Polygon(
            points={{-60,10},{-60,20},{-40,40},{-40,-40},{-60,-20},{-60,10}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={128,128,128}),
          Polygon(
            points={{60,20},{40,40},{40,-40},{60,-20},{60,20}},
            lineColor={128,128,128},
            fillColor={128,128,128},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,-90},{-50,-90},{-20,-30},{20,-30},{48,-90},{60,-90},{60,
                -100},{-60,-100},{-60,-90}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Text(
          extent={{-150,110},{150,70}},
          textString="%name",
          lineColor={0,0,255},
          fontSize=0),
          Text(
          extent={{-150,-110},{150,-150}},
          lineColor={0,0,0},
          textString="ratio=%ratio")}),
      Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics),
      conversion(noneFromVersion=""));
  end Gear;

  model ControlledMotor "Current controlled DC motor"
    extends Modelica.Electrical.Machines.Icons.Machine;

    parameter Real km=30 "Gain";
    parameter Modelica.SIunits.Time Tm=0.005 "Time Constant (T>0 required)";

    Modelica.Blocks.Continuous.PI PI(k=km, T=Tm)
      annotation (Placement(transformation(extent={{-52,60},{-32,80}},
          rotation=0)));

    Modelica.Mechanics.Rotational.Interfaces.Flange_b flange
      annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=
           0)));
    Modelica.Blocks.Interfaces.RealInput refCurrent
    "Reference current of motor"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0)));
    Modelica.Blocks.Math.Feedback feedback
      annotation (Placement(transformation(extent={{-90,-10},{-70,10}},
          rotation=0)));
    Modelica.Blocks.Continuous.FirstOrder firstOrder(T=0.001, initType=Modelica.Blocks.Types.Init.InitialState)
      annotation (Placement(transformation(extent={{-18,60},{2,80}},  rotation=
            0)));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=13.8)
      annotation (Placement(transformation(extent={{-10,20},{10,40}},  rotation=0)));
    Modelica.Electrical.Analog.Basic.Inductor inductor(L=0.061, i(start=0,
        fixed=true))
      annotation (Placement(transformation(extent={{20,20},{40,40}},rotation=0)));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{40,-50},{60,-30}}, rotation=0)));
    Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage
      annotation (Placement(transformation(
          origin={-20,0},
          extent={{-10,10},{10,-10}},
          rotation=270)));
    Modelica.Electrical.Analog.Basic.EMF emf(k=1.016)
                                             annotation (Placement(transformation(
            extent={{40,-10},{60,10}},
                                     rotation=0)));
    Modelica.Mechanics.Rotational.Components.Inertia motorInertia(J=0.0025)
      annotation (Placement(transformation(extent={{66,-10},{86,10}},
                                                                    rotation=0)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
      annotation (Placement(transformation(extent={{18,-40},{-2,-20}},  rotation=
              0)));
  equation

    connect(feedback.u1, refCurrent)
      annotation (Line(points={{-88,0},{-120,0}}, color={0,0,127}));
    connect(feedback.y, PI.u)
      annotation (Line(points={{-71,0},{-66,0},{-66,70},{-54,70}},
                                                 color={0,0,127}));
    connect(PI.y, firstOrder.u)
      annotation (Line(points={{-31,70},{-20,70}},
                                                 color={0,0,127}));
    connect(signalVoltage.p,resistor. p) annotation (Line(points={{-20,10},{-20,
          30},{-10,30}},   color={0,0,255}));
    connect(resistor.n,inductor. p)
      annotation (Line(points={{10,30},{20,30}}, color={0,0,255}));
    connect(inductor.n,emf. p) annotation (Line(points={{40,30},{50,30},{50,10}},
          color={0,0,255}));
    connect(emf.n,ground. p)
      annotation (Line(points={{50,-10},{50,-30}},
                                                 color={0,0,255}));
    connect(motorInertia.flange_a,emf.flange)
      annotation (Line(points={{66,0},{60,0}},   color={0,0,0}));
    connect(currentSensor.n,signalVoltage. n) annotation (Line(points={{-2,-30},
          {-20,-30},{-20,-10}}, color={0,0,255}));
    connect(currentSensor.p,ground. p)
      annotation (Line(points={{18,-30},{50,-30}}, color={0,0,255}));
  connect(firstOrder.y, signalVoltage.v) annotation (Line(points={{3,70},{28,70},
          {28,46},{-50,46},{-50,0},{-32,0}}, color={0,0,127}));
  connect(currentSensor.i, feedback.u2) annotation (Line(points={{8,-41},{8,-52},
          {-80,-52},{-80,-8}}, color={0,0,127}));
  connect(motorInertia.flange_b, flange)
    annotation (Line(points={{86,0},{100,0}}, color={0,0,0}));
    annotation (         Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
        Text(
          extent={{-150,120},{150,80}},
          lineColor={0,0,255},
          textString="%name"),
          Rectangle(
            extent={{-100,10},{-60,-10}},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={95,95,95})}));
  end ControlledMotor;

  model SpeedController "PI Geschwindigkeits-Regler"
    import SI = Modelica.SIunits;
    parameter Real ks "Verstaerkung vom PI Geschwindigkeitsregler";
    parameter SI.Time Ts "Zeitkonstante vom PI Geschwindigkeitsregler";
    parameter Real ratio=105 "Getriebe-Uebersetzung";

    Modelica.Blocks.Interfaces.RealInput refSpeed "Reference speed of load"
      annotation ( Placement(
        transformation(extent={{-140,-20},{-100,20}}, rotation=0)));
    Modelica.Blocks.Interfaces.RealOutput refCurrent "Reference current of motor"
      annotation (Placement(
        transformation(extent={{100,-10},{120,10}}, rotation=0)));
    Modelica.Blocks.Interfaces.RealInput motorSpeed "Actual speed of motor"
                                                      annotation (
      Placement(transformation(
        origin={0,-120},
        extent={{-20,-20},{20,20}},
        rotation=90)));
    Modelica.Blocks.Math.Feedback feedback
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, rotation=
           0)));
    Modelica.Blocks.Math.Gain gain(k=ratio)
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}},
          rotation=0)));
    Modelica.Blocks.Continuous.PI PI(T=Ts, k=ks)
      annotation (Placement(transformation(extent={{30,-10},{50,10}}, rotation=
            0)));
  equation
    connect(gain.y,feedback.u1)             annotation (Line(points={{-29,0},{
          -8,0}}, color={0,0,127}));
    connect(feedback.y,PI.u)             annotation (Line(points={{9,0},{28,0}},
        color={0,0,127}));
    connect(PI.y,       refCurrent) annotation (Line(points={{51,0},{110,0}},
        color={0,0,127}));
    connect(gain.u,      refSpeed)
      annotation (Line(points={{-52,0},{-120,0}}, color={0,0,127}));
    connect(feedback.u2,      motorSpeed)
      annotation (Line(points={{0,-8},{0,-120}}, color={0,0,127}));
    annotation (
      Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
      Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
          Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,0},
            pattern=LinePattern.Solid,
            lineThickness=0.25,
            fillColor={235,235,235},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-30,54},{30,24}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-30,40},{-60,50},{-60,30},{-30,40}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-31,-41},{-78,-41},{-78,39},{-30,39}}),
          Rectangle(
            extent={{-30,-26},{30,-56}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{60,-32},{30,-42},{60,-52},{60,-32}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Line(points={{30,39},{76,39},{76,-41},{30,-41}}),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={0,0,255},
            fontSize=0)}),
      conversion(noneFromVersion=""));
  end SpeedController;

  model Servo "Drehzahlgeregelter Motor mit Getriebe"
    extends Modelica.Electrical.Machines.Icons.TransientMachine;
    Modelica.Blocks.Interfaces.RealInput refSpeed "Reference speed of load" annotation (Placement(
        transformation(extent={{-140,-20},{-100,20}}, rotation=0)));
    Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b annotation (Placement(
        transformation(extent={{90,-10},{110,10}}, rotation=0)));

    parameter Real ks "Verstaerkung vom PI Geschwindigkeitsregler";
    parameter Modelica.SIunits.Time Ts "Zeitkonstante vom PI Geschwindigkeitsregler";

    parameter Real km=30 "Verstaerkung vom PI Motorregler";
    parameter Modelica.SIunits.Time Tm=0.005 "Zeitkonstante vom PI Motorregler";
    parameter Real ratio=105 "Getriebeuebersetzung";

    ServoSystem.ControlledMotor motor(km=km, Tm=Tm) annotation (Placement(
          transformation(extent={{-20,-10},{0,10}}, rotation=0)));
    SpeedController speedController(
      ks=ks,
      Ts=Ts,
      ratio=ratio) annotation (Placement(transformation(extent={{-70,-10},{-50,10}},
            rotation=0)));
    Gear gear(ratio=ratio) annotation (Placement(transformation(extent={{30,-10},
              {50,10}}, rotation=0)));

    Modelica.Blocks.Math.Feedback speedError
      annotation (Placement(transformation(extent={{-40,40},{-20,60}}, rotation=
           0)));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor2
      annotation (Placement(transformation(
        origin={70,20},
        extent={{-10,-10},{10,10}},
        rotation=90)));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor1
      annotation (Placement(transformation(
        origin={10,-20},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  equation
    connect(speedSensor2.w, speedError.u2)
      annotation (Line(points={{70,31},{70,34},{-30,34},{-30,42}}, color={0,0,
          127}));
    connect(refSpeed, speedController.refSpeed)
      annotation (Line(points={{-120,0},{-72,0}}, color={0,0,127}));
    connect(speedController.refCurrent, motor.refCurrent)
      annotation (Line(points={{-49,0},{-22,0}}, color={0,0,127}));
    connect(speedError.u1,      refSpeed)
      annotation (Line(points={{-38,50},{-90,50},{-90,0},{-120,0}}, color={0,0,
          127}));
    connect(motor.flange,speedSensor1.flange)    annotation (Line(points={{0,0},
          {10,0},{10,-10}}, color={0,0,0}));
    connect(speedSensor1.w, speedController.motorSpeed) annotation (Line(points=
           {{10,-31},{10,-40},{-60,-40},{-60,-12}}, color={0,0,127}));
    connect(motor.flange, gear.flange_a)
      annotation (Line(points={{0,0},{30,0}}, color={0,0,0}));
    connect(gear.flange_b, flange_b)
      annotation (Line(points={{50,0},{100,0}}, color={0,0,0}));
    connect(gear.flange_b, speedSensor2.flange)
      annotation (Line(points={{50,0},{70,0},{70,10}}, color={0,0,0}));
    annotation (
      Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics={
        Text(
          extent={{-150,120},{150,80}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-150,-110},{150,-150}},
          lineColor={0,0,0},
          textString="ks=%ks, Ts=%Ts"),
          Rectangle(
            extent={{-100,10},{-60,-10}},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={95,95,95})}),
      Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics));
  end Servo;

model Test "Test to tune speed controller"
  import SI = Modelica.SIunits;
  extends Modelica.Icons.Example;
  parameter Real ks = 0.8 "Verstaerkung vom PI Geschwindigkeitsregler";
  parameter SI.Time Ts= 0.08 "Zeitkonstante vom PI Geschwindigkeitsregler";
    Servo servo(ks=ks, Ts=Ts) annotation (Placement(transformation(extent={{-20,
              0},{0,20}}, rotation=0)));

  Modelica.Mechanics.Rotational.Components.Inertia load(J=170)
      "J=50..170 kgm^2" annotation (Placement(transformation(extent={{20,0},{40,
              20}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp(duration=1.18, height=2.95)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}}, rotation=
          0)));
equation
    connect(ramp.y, servo.refSpeed)
      annotation (Line(points={{-39,10},{-22,10}}, color={0,0,127}));
    connect(servo.flange_b, load.flange_a)
      annotation (Line(points={{0,10},{20,10}}, color={0,0,0}));
  annotation (
    Diagram(coordinateSystem(
      preserveAspectRatio=false,
      extent={{-100,-100},{100,100}},
      grid={2,2}), graphics),
    experiment(StopTime=2),
  __Dymola_Commands(file="Scripts/plot load.w, speedError1.mos"
      "plot load.w, speedError", file="Scripts/plot current1.mos"
      "plot current"));
end Test;
  annotation (uses(Modelica(version="3.2.3")));
end ServoSystem;
