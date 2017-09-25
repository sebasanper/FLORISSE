# -*- coding: utf-8 -*-


class neutral:
    """A control set for the turbines in the windfarm"""
    def __init__(self, layout):
        nTurbs = layout.nTurbs

        # individual turbine parameters (update Ct and Cp)
        self.yawAngles = [0 for i in range(nTurbs)]     # yaw angles [deg]
        self.tiltAngles = [0 for i in range(nTurbs)]    # tilt angles [deg]
        self.bladePitch = [1.9 if turb.usePitch else 0
                           for turb in layout.turbines]  # blade pitch [deg]
        self.TSR = [8.0 if turb.useTSR else 0
                    for turb in layout.turbines]  # Tip speed ratio [-]


class yawed(neutral):
    def __init__(self, layout):
        super().__init__(layout)
        self.yawAngles[0] = 20
