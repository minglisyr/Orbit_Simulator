def main():
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    import bpy

    __author__ = 'LI, MING'
    z = int(864000)
    G = 6.67e-11
    angleCF = math.pi / 180
    pi = math.pi
    twopi = 2 * math.pi

    t = [None]*z
    t[0] = 0
    dt = 1

    class Mass:
        def __init__(self, mass, x0 ,y0):
            self.mass = mass
            self.mu = G * self.mass
            self.x = [None]*z
            self.y = [None]*z
            self.r = [None]*z
            self.x[0] = x0
            self.y[0] = y0
            self.r[0] = math.sqrt(x0**2 + y0**2)

    class CelestialBody(Mass):
        def Settings(self, radius, period, dE):
            self.radius = radius
            self.period = period * 24 * 60 * 60.0
            self.dE = dE
            self.soi = (self.mass / Earth.mass) ** 0.4 * self.radius
        def Position(self):
            self.x[i+1] = (self.dE) * \
            math.cos(2 * math.pi / self.period * t[i+1] + self.omega[0])
            self.y[i+1] = (self.dE) * \
            math.sin(2 * math.pi / self.period * t[i+1] + self.omega[0])
            self.r[i+1] = math.sqrt(self.x[i+1]**2 + self.y[i+1]**2)
        def MoonAngle_init(self):
            self.omega = [None]*z
            if self.y[0] >= 0:
                self.omega[0] = math.acos(self.x[0] / self.r[0])
            elif self.y[0] < 0:
                self.omega[0] = twopi - math.acos(self.x[0] / self.r[0])
            else:
                print("omega error in CelestialBody.MoonAngle_init()")
        def Angles(self):
            if self.y[i+1] >= 0:
                self.omega[i+1] = math.acos(self.x[i+1] / self.r[i+1])
            elif self.y[i+1] < 0:
                self.omega[i+1] = twopi - math.acos(self.x[i+1] / self.r[i+1])
            else:
                print("omega error in CelestialBody.Angles()")
        
    Earth = CelestialBody(5.97237e+24, 0.0, 0.0)
    Earth.Settings(6.3781e+6, 365, 0.0)
    Moon = CelestialBody(7.342e+22, 364658999.01580966, 121599365.28529961)
    Moon.Settings(1.7371e+6, 30, 3.84399e+8)
    Moon.MoonAngle_init()

    class Rocket(Mass):
        def Settings(self, crash, procedureTurnTime, procedureTurnAngle, \
                    Turnoff, TargetAltitude, TargetAltitude2, \
                    TargetMoonAltitude, v_x0, v_y0, theta, epsilon, Fr):
            self.crash = crash
            self.procedureTurnTime = procedureTurnTime
            self.procedureTurnAngle = procedureTurnAngle * angleCF
            self.Turnoff = Turnoff
            self.TargetAltitude = TargetAltitude
            self.TargetAltitude2 = TargetAltitude2
            self.TargetMoonAltitude = TargetMoonAltitude
            self.Alt_Earth = [None]*z
            self.Mx = [None]*z
            self.My = [None]*z
            self.Mr = [None]*z
            self.Alt_Moon = [None]*z
            self.v_x = [None]*z
            self.v_y = [None]*z
            self.v = [None]*z
            self.theta = [None]*z
            self.epsilon = [None]*z
            self.tau = [None]*z
            self.phi = [None]*z
            self.Fg_Earth = [None]*z
            self.Fg_Moon = [None]*z
            self.Fr = [None]*z
            self.a_x = [None]*z
            self.a_y = [None]*z
            self.a = [None]*z
            self.Alt_Earth[0] = 0
            self.Mx[0] = self.x[0] - Moon.x[0]
            self.My[0] = self.y[0] - Moon.y[0]
            self.Mr[0] = math.sqrt(self.Mx[0]**2 + self.My[0]**2)
            self.Alt_Moon[0] = self.Mr[0] - Moon.radius
            self.v_x[0] = v_x0
            self.v_y[0] = v_y0
            self.v[0] = math.sqrt(self.v_x[0]**2 + self.v_y[0]**2)
            self.theta[0] = theta
            self.epsilon[0] = epsilon
            self.tau[0] = self.epsilon[0]
            if self.y[0] >= Moon.y[0]:
                self.phi[0] = math.acos(self.Mx[0] / self.Mr[0])
            elif self.y[0] < Moon.y[0]:
                self.phi[0] = twopi - math.acos(self.Mx[0] / self.Mr[0])
            else:
                print("phi error in Rocket.Settings()", t[0])
            self.Fg_Earth[0] = -G * Earth.mass * self.mass / self.r[0]**2
            self.Fg_Moon[0] = -G * Moon.mass * self.mass / self.Mr[0]**2
            self.Fr[0] = Fr
            self.a_x[0] = (self.Fg_Earth[0] * math.cos(self.theta[0]) + \
                    self.Fr[0] * math.cos(self.epsilon[0])) / self.mass
            self.a_y[0] = (self.Fg_Earth[0] * math.sin(self.theta[0]) + \
                    self.Fr[0] * math.sin(self.epsilon[0])) / self.mass
            self.a[0] = math.sqrt(self.a_x[0]**2 + self.a_y[0]**2)
            self.Acceleration_ProcedureTurn_called = 0
            self.Thrust_Earth_InCircle_called = 0
            self.Thrust_Earth_InEllipse_called = 0
            self.Thrust_Moon_InCircle_called = 0
        def Position(self):
            self.x[i+1] = self.x[i] + self.v_x[i] * dt + 0.5 * self.a_x[i] * dt**2
            self.y[i+1] = self.y[i] + self.v_y[i] * dt + 0.5 * self.a_y[i] * dt**2
            self.r[i+1] = math.sqrt(self.x[i+1]**2 + self.y[i+1]**2)
            self.Alt_Earth[i+1] = self.r[i+1] - Earth.radius
            self.Mx[i+1] = self.x[i+1] - Moon.x[i+1]
            self.My[i+1] = self.y[i+1] - Moon.y[i+1]
            self.Mr[i+1] = math.sqrt(self.Mx[i+1]**2 + self.My[i+1]**2)
            self.Alt_Moon[i+1] = self.Mr[i+1] - Moon.radius
            while i in keyframe:
                temp_Moon_x = Moon.x[i+1] / Earth.radius
                temp_Moon_y = Moon.y[i+1] / Earth.radius
                temp_Stamina_x = Stamina.x[i+1] / Earth.radius
                temp_Stamina_y = Stamina.y[i+1] / Earth.radius
                temp_Moon_x_list.append(temp_Moon_x)
                temp_Moon_y_list.append(temp_Moon_y)
                temp_Stamina_x_list.append(temp_Stamina_x)
                temp_Stamina_y_list.append(temp_Stamina_y)
                break
        def Velocity(self):
            self.v_x[i+1] = self.v_x[i] + self.a_x[i] * dt
            self.v_y[i+1] = self.v_y[i] + self.a_y[i] * dt
            self.v[i+1] = math.sqrt(self.v_x[i+1]**2 + self.v_y[i+1]**2)
        def Force(self):
            self.Fg_Earth[i+1] = -G * Earth.mass * self.mass / self.r[i+1]**2
            self.Fg_Moon[i+1] = -G * Moon.mass * self.mass / self.Mr[i+1]**2
            if t[i+1] < Stamina.Turnoff:
                self.Fr[i+1] = self.Fr[0]
            else:
                self.Fr[i+1] = 0
        def Angles(self):   
            if self.v_y[i+1] >= 0:
                self.epsilon[i+1] = math.acos(self.v_x[i+1] / self.v[i+1])
            elif self.v_y[i+1] < 0:
                self.epsilon[i+1] = twopi - math.acos(self.v_x[i+1] / self.v[i+1])
            else:
                print("epsilon error in Rocket.Angles()", t[i+1])
            if self.y[i+1] >= 0:
                self.theta[i+1] = math.acos(self.x[i+1] / self.r[i+1])
            elif self.y[i+1] < 0:
                self.theta[i+1] = twopi - math.acos(self.x[i+1] / self.r[i+1])
            else:
                print("theta error in Rocket.Angles()", t[i+1])
            if self.y[i+1] >= Moon.y[i+1]:
                self.phi[i+1] = math.acos(self.Mx[i+1] / self.Mr[i+1])
            elif self.y[i+1] < Moon.y[i+1]:
                self.phi[i+1] = twopi - math.acos(self.Mx[i+1] / self.Mr[i+1])
            else:
                print("phi error in Rocket.Angles()", t[i+1])
        def Acceleration_Normal(self):
            self.tau[i+1] = self.epsilon[i+1]
            self.a_x[i+1] = (self.Fg_Earth[i+1] * math.cos(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.cos(self.phi[i+1]) + \
                    self.Fr[i+1] * math.cos(self.tau[i+1])) / self.mass
            self.a_y[i+1] = (self.Fg_Earth[i+1] * math.sin(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.sin(self.phi[i+1]) + \
                    self.Fr[i+1] * math.sin(self.tau[i+1])) / self.mass
            self.a[i+1] = math.sqrt(self.a_x[i+1]**2 + self.a_y[i+1]**2)
        def Acceleration_ProcedureTurn(self):
            self.Acceleration_ProcedureTurn_called += 1
            a_tx = self.v[i+1] * math.cos(self.procedureTurnAngle) - self.v_x[i+1]
            a_ty = self.v[i+1] * math.sin(self.procedureTurnAngle) - self.v_y[i+1]
            a_t = math.sqrt(a_tx**2 + a_ty**2)
            if a_ty >= 0: 
                self.tau[i+1] = math.acos(a_tx / a_t)
            elif a_ty < 0:
                self.tau[i+1] = twopi - math.acos(a_tx / a_t)
            else:
                print("tau error in Rocket.Acceleration_ProcedureTurn()", t[i+1])
            self.a_x[i+1] = (self.Fg_Earth[i+1] * math.cos(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.cos(self.phi[i+1]) + \
                    self.Fr[i+1] * math.cos(self.tau[i+1])) / self.mass
            self.a_y[i+1] = (self.Fg_Earth[i+1] * math.sin(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.sin(self.phi[i+1]) + \
                    self.Fr[i+1] * math.sin(self.tau[i+1])) / self.mass
            self.a[i+1] = math.sqrt(self.a_x[i+1]**2 + self.a_y[i+1]**2)
        def Thrust_Earth_InCircle(self):
            self.Thrust_Earth_InCircle_called += 1
            gamma = pi / 2 + self.theta[i+1]
            v_o = math.sqrt(Earth.mu * self.r[i+1]) / self.r[i+1]
            a_tx = (v_o * math.cos(gamma) - self.v_x[i+1]) / dt
            a_ty = (v_o * math.sin(gamma) - self.v_y[i+1]) / dt
            a_t = math.sqrt(a_tx**2 + a_ty**2)
            if a_ty >= 0:
                self.tau[i+1] = math.acos(a_tx / a_t)
            elif a_ty < 0:
                self.tau[i+1] = twopi - math.acos(a_tx / a_t)
            else:
                print("tau error in Rocket.Thrust_Earth_InCircle()", t[i+1])
            self.Fr[i+1] = self.mass * a_t
            self.a_x[i+1] = (self.Fg_Earth[i+1] * math.cos(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.cos(self.phi[i+1]) + \
                    self.Fr[i+1] * math.cos(self.tau[i+1])) / self.mass
            self.a_y[i+1] = (self.Fg_Earth[i+1] * math.sin(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.sin(self.phi[i+1]) + \
                    self.Fr[i+1] * math.sin(self.tau[i+1])) / self.mass
            self.a[i+1] = math.sqrt(self.a_x[i+1]**2 + self.a_y[i+1]**2)
            print('In circular orbit at', self.Thrust_Earth_InCircle_called, t[i+1], 's')
        def Thrust_Earth_InEllipse(self, r_a):
            self.Thrust_Earth_InEllipse_called += 1
            gamma = pi / 2 + self.theta[i+1]
            v_o = math.sqrt(2 * Earth.mu * \
                (r_a + Earth.radius) * (self.r[i+1]) / \
                ((r_a + Earth.radius) + (self.r[i+1]))) / \
                self.r[i+1]
            a_tx = (v_o * math.cos(gamma) - self.v_x[i+1]) / dt
            a_ty = (v_o * math.sin(gamma) - self.v_y[i+1]) / dt
            a_t = math.sqrt(a_tx**2 + a_ty**2)
            if a_ty >= 0:
                self.tau[i+1] = math.acos(a_tx / a_t)
            elif a_ty < 0:
                self.tau[i+1] = twopi - math.acos(a_tx / a_t)
            else:
                print("tau error in Rocket.Thrust_Earth_InEllipse()", t[i+1])
            self.Fr[i+1] = self.mass * a_t
            self.a_x[i+1] = (self.Fg_Earth[i+1] * math.cos(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.cos(self.phi[i+1]) + \
                    self.Fr[i+1] * math.cos(self.tau[i+1])) / self.mass
            self.a_y[i+1] = (self.Fg_Earth[i+1] * math.sin(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.sin(self.phi[i+1]) + \
                    self.Fr[i+1] * math.sin(self.tau[i+1])) / self.mass
            self.a[i+1] = math.sqrt(self.a_x[i+1]**2 + self.a_y[i+1]**2)
            print('In transfer orbit', self.Thrust_Earth_InEllipse_called, t[i+1], 's')
        def Thrust_Moon_InCircle(self):
            self.Thrust_Moon_InCircle_called += 1
            gamma = self.phi[i+1] - pi / 2
            Lambda = Moon.omega[i+1] + pi / 2
            v_m = math.sqrt(Earth.mass * G / Moon.dE)
            v_o = math.sqrt(Moon.mu * self.Mr[i+1]) / self.Mr[i+1]
            a_tx = (v_o * math.cos(gamma) + v_m * \
                math.cos(Lambda) - self.v_x[i+1]) / dt
            a_ty = (v_o * math.sin(gamma) + v_m * \
                math.sin(Lambda) - self.v_y[i+1]) / dt
            a_t = math.sqrt(a_tx**2 + a_ty**2)
            if a_ty >= 0:
                self.tau[i+1] = math.acos(a_tx / a_t)
            elif a_ty < 0:
                self.tau[i+1] = twopi - math.acos(a_tx / a_t)
            else:
                print("tau error in Rocket.Thrust_Moon_InCircle()", t[i+1])
            self.Fr[i+1] = self.mass * a_t
            self.a_x[i+1] = (self.Fg_Earth[i+1] * math.cos(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.cos(self.phi[i+1]) + \
                    self.Fr[i+1] * math.cos(self.tau[i+1])) / self.mass
            self.a_y[i+1] = (self.Fg_Earth[i+1] * math.sin(self.theta[i+1]) + \
                    self.Fg_Moon[i+1] * math.sin(self.phi[i+1]) + \
                    self.Fr[i+1] * math.sin(self.tau[i+1])) / self.mass
            self.a[i+1] = math.sqrt(self.a_x[i+1]**2 + self.a_y[i+1]**2)
        def Acceleration(self):
            if t[i+1] > self.procedureTurnTime and \
            self.procedureTurnAngle - self.epsilon[i] > 0 and \
            self.Thrust_Earth_InCircle_called == 0:
                self.Acceleration_ProcedureTurn()
            elif t[i+1] > self.procedureTurnTime and \
            self.TargetAltitude + 5e5 > self.Alt_Earth[i+1] > \
            self.TargetAltitude and \
            self.Thrust_Earth_InCircle_called == 0 and \
            self.Thrust_Earth_InEllipse_called == 0:
                self.Thrust_Earth_InCircle()
            elif t[i+1] > self.procedureTurnTime and \
            self.TargetAltitude + 5e5 > self.Alt_Earth[i+1] > \
            self.TargetAltitude - 5e5 and \
            self.x[i+1] >= 0 and self.y[i+1] < 0 and \
            self.Thrust_Earth_InEllipse_called == 0 and \
            self.Thrust_Earth_InCircle_called == 1:
                self.Thrust_Earth_InEllipse(self.TargetAltitude2)
            elif t[i+1] > self.procedureTurnTime and \
            self.TargetAltitude2 + 5e5 > self.Alt_Earth[i+1] > \
            self.TargetAltitude2 and \
            self.Thrust_Earth_InCircle_called == 1 and \
            self.Thrust_Earth_InEllipse_called == 1:
                self.Thrust_Earth_InCircle()
            elif t[i+1] > self.procedureTurnTime and \
            self.TargetAltitude2 + 5e5 > self.Alt_Earth[i+1] > \
            self.TargetAltitude2 - 5e5 and \
            self.x[i+1] >= 0 and self.y[i+1] < 0 and \
            self.Thrust_Earth_InCircle_called == 2 and \
            self.Thrust_Earth_InEllipse_called == 1:
                self.Thrust_Earth_InEllipse(Moon.dE-Earth.radius)
            elif t[i+1] > self.procedureTurnTime and \
            225 * angleCF > self.epsilon[i] > 180 * angleCF and \
            Moon.dE + 5e5 > self.Alt_Earth[i+1] > Moon.dE and \
            self.Thrust_Earth_InCircle_called == 2 and \
            self.Thrust_Moon_InCircle_called == 0:
                self.Thrust_Earth_InCircle()
            elif t[i+1] > self.procedureTurnTime and \
            self.TargetMoonAltitude > self.Alt_Moon[i+1]:
                self.Thrust_Moon_InCircle()
            else:
                self.Acceleration_Normal() 
        def Crash(self):
            if self.Alt_Earth[i+1] <= 0.0 or self.Alt_Moon[i+1] <= 0.0:
                print("crash at", t[i+1], "s")
                self.crash = True
            else:
                pass
    Stamina = Rocket(1.0e+5, Earth.radius, 0.0)
    Stamina.Settings(False, 20, 45, 200, 1e6, 10e6, 1e6, 0.0, 0.0, 0.0, 0.001*angleCF, 5.5e6)
    temp_Moon_x_list = []
    temp_Moon_y_list = []
    temp_Stamina_x_list = []
    temp_Stamina_y_list = []
    temp_Moon_x_list.append(Moon.x[0] / Earth.radius)
    temp_Moon_y_list.append(Moon.y[0] / Earth.radius)
    temp_Stamina_x_list.append(Stamina.x[0] / Earth.radius)
    temp_Stamina_y_list.append(Stamina.y[0] / Earth.radius)

    keyframe = np.arange(0, z, 500)

    for i in range(z-1):
        t[i+1] = t[i] + dt
        Moon.Position()
        Moon.Angles()
        Stamina.Position()
        Stamina.Velocity()
        Stamina.Force()
        Stamina.Angles()
        Stamina.Acceleration()
        Stamina.Crash()
        if Stamina.crash:
            break
        else:
            pass
    
    open("blender_data_list.txt", \
     "w").write('temp_Moon_x_list = ' + str(temp_Moon_x_list) + '\n' \
     + 'temp_Moon_y_list = ' + str(temp_Moon_y_list) + '\n' \
     + 'temp_Stamina_x_list = ' + str(temp_Stamina_x_list) + '\n' \
     + 'temp_Stamina_y_list = ' +  str(temp_Stamina_y_list) + '\n')    
    
    bpy.data.scenes["Scene"].frame_current = 0
    for n in range(int(z/500)):
        bpy.ops.anim.keyframe_insert_menu(type='Location')
        bpy.data.objects["Moon"].location[0] = temp_Moon_x_list[n]
        bpy.data.objects["Moon"].location[1] = temp_Moon_y_list[n]
        bpy.data.objects["CameraMoon"].location[0] = temp_Moon_x_list[n]
        bpy.data.objects["CameraMoon"].location[1] = temp_Moon_y_list[n]
        bpy.data.objects["CameraStamina"].location[0] = temp_Stamina_x_list[n]
        bpy.data.objects["CameraStamina"].location[1] = temp_Stamina_y_list[n]
        bpy.data.objects["Stamina"].location[0] = temp_Stamina_x_list[n]
        bpy.data.objects["Stamina"].location[1] = temp_Stamina_y_list[n]
        bpy.data.scenes["Scene"].frame_current += 1
    
    def graph_Common_Settings():
        plt.xlim(xmin=0)
        plt.xlim(xmax=t[i+1])
        plt.hlines(0, t[0], t[i+1], "gray", linewidth=0.5)

    def graph_Position():
        plt.figure(figsize=(8,5), dpi = 300)
        plt.subplot(221)
        plt.title('t, Stamina.x')
        plt.plot(t, Stamina.x, "orange")
        graph_Common_Settings()
        plt.subplot(222)
        plt.title('t, Stamina.y')
        plt.plot(t, Stamina.y, "orange")
        graph_Common_Settings()
        plt.subplot(223)
        plt.title('t, Stamina.Alt_Earth')
        plt.plot(t, Stamina.Alt_Earth, "orange")
        plt.hlines(Stamina.TargetAltitude, t[0], t[i+1], "gray", linewidth=0.5)
        plt.hlines(Stamina.TargetAltitude2, t[0], t[i+1], "gray", linewidth=0.5)
        graph_Common_Settings()
        plt.subplot(224)
        plt.title('t, Stamina.Alt_Moon')
        plt.plot(t, Stamina.Alt_Moon, "orange")
        plt.hlines(Stamina.TargetMoonAltitude, t[0], t[i+1], "gray", linewidth=0.5)
        graph_Common_Settings()
        plt.tight_layout(w_pad=2.0, h_pad=2.0)
        plt.show()
    graph_Position()

    def graph_Velocity():
        plt.figure(figsize=(8,5), dpi = 300)
        plt.subplot(221)
        plt.title('t, Stamina.v_x')
        plt.plot(t, Stamina.v_x, "green")
        graph_Common_Settings()
        plt.subplot(222)
        plt.title('t, Stamina.v_y')
        plt.plot(t, Stamina.v_y, "green")
        graph_Common_Settings()
        plt.subplot(223)
        plt.title('t, Stamina.v')
        plt.plot(t, Stamina.v, "green")
        graph_Common_Settings()
        plt.tight_layout(w_pad=2.0, h_pad=2.0)
        plt.show()
    graph_Velocity()

    def graph_Acceleration():
        plt.figure(figsize=(8,5), dpi = 300)
        plt.subplot(221)
        plt.title('t, Stamina.a_x')
        plt.plot(t, Stamina.a_x, "red")
        graph_Common_Settings()
        plt.subplot(222)
        plt.title('t, Stamina.a_y')
        plt.plot(t, Stamina.a_y, "red")
        graph_Common_Settings()
        plt.subplot(223)
        plt.title('t, Stamina.a')
        plt.plot(t, Stamina.a, "red")
        graph_Common_Settings()
        plt.tight_layout(w_pad=2.0, h_pad=2.0)
        plt.show()
    graph_Acceleration()

    def graph_Force():
        plt.figure(figsize=(8,5), dpi = 300)
        plt.subplot(221)
        plt.title('t, Stamina.Fg_Earth')
        plt.plot(t, Stamina.Fg_Earth, "purple")
        graph_Common_Settings()
        plt.subplot(222)
        plt.title('t, Stamina.Fg_Moon')
        plt.plot(t, Stamina.Fg_Moon, "purple")
        graph_Common_Settings()
        plt.subplot(223)
        plt.title('t, Stamina.Fr')
        plt.plot(t, Stamina.Fr, "purple")
        graph_Common_Settings()
        plt.tight_layout(w_pad=2.0, h_pad=2.0)
        plt.show()
    graph_Force()

    def graph_Angles_Settings():
        plt.hlines(0, t[0], t[i+1], "gray", linewidth=0.5)
        plt.hlines(pi/2, t[0], t[i+1], "gray", linewidth=0.5)
        plt.hlines(pi, t[0], t[i+1], "gray", linewidth=0.5)
        plt.hlines(pi * 3 / 2, t[0], t[i+1], "gray", linewidth=0.5)
        plt.hlines(twopi, t[0], t[i+1], "gray", linewidth=0.5)
        plt.xlim(xmin=0)
        plt.xlim(xmax=t[i+1])
        plt.ylim(ymin=0)
        plt.ylim(ymax=twopi)

    def graph_Angles():
        plt.figure(figsize=(8,5), dpi = 300)
        plt.subplot(221)
        plt.title('t, Stamina.epsilon')
        plt.plot(t, Stamina.epsilon, ',', color="blue")
        graph_Angles_Settings()
        plt.subplot(222)
        plt.title('t, Stamina.theta')
        plt.plot(t, Stamina.theta, ',', color="blue")
        graph_Angles_Settings()
        plt.subplot(223)
        plt.title('t, Stamina.phi')
        plt.plot(t, Stamina.phi, ',', color="blue")
        graph_Angles_Settings()
        plt.subplot(224)
        plt.title('t, Stamina.tau')
        plt.plot(t, Stamina.tau, ',', color="blue")
        graph_Angles_Settings()
        plt.tight_layout(w_pad=2.0, h_pad=2.0)
        plt.show()
    graph_Angles()

    def graph_EarthProximity():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Earth in Proximity')
        plt.plot(Stamina.x, Stamina.y, color="orange")
        #plt.plot(Interstellar.x, Interstellar.y, color="brown")
        plt.xlim(xmin=-Earth.radius-3e7)
        plt.xlim(xmax=Earth.radius+3e7)
        plt.ylim(ymin=-Earth.radius-3e7)
        plt.ylim(ymax=Earth.radius+3e7)
        plt.hlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        plt.vlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((0,0), Stamina.TargetAltitude + Earth.radius, 
                                color="blue", fill=None))
        ax.add_patch(plt.Circle((0,0), Stamina.TargetAltitude2 + Earth.radius, 
                                color="blue", fill=None))
        ax.add_patch(plt.Circle((Stamina.x[t.index(max(t))],
                                        Stamina.y[t.index(max(t))]), 
                                5e4, color="orange", fill=True))
        plt.show()
    graph_EarthProximity()

    def graph_MoonProximity():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Moon in Proximity')
        plt.plot(Stamina.x, Stamina.y, color="orange")
        plt.plot(Moon.x, Moon.y, color="#377EB8")
        plt.xlim(xmin=Moon.x[i+1]-5.0e6)
        plt.xlim(xmax=Moon.x[i+1]+5.0e6)
        plt.ylim(ymin=Moon.y[i+1]-5.0e6)
        plt.ylim(ymax=Moon.y[i+1]+5.0e6)
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((Stamina.x[t.index(max(t))],
                                        Stamina.y[t.index(max(t))]),
                                5e4, color="orange", fill=True))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), 
                                5e4, color="#377EB8", fill=True))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), 
                                Moon.radius, color="#377EB8", fill=None))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), 
                                Moon.radius + Stamina.TargetMoonAltitude, 
                                color="purple", fill=None))
        plt.show()
    graph_MoonProximity()

    def graph_EarthSurface():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Earth Surface')
        plt.plot(Stamina.x, Stamina.y, color="orange")
        plt.xlim(xmin=Earth.radius-1e4)
        plt.xlim(xmax=Earth.radius+1e6*2)
        plt.ylim(ymin=-1e6)
        plt.ylim(ymax=1e6)
        plt.hlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        plt.vlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((Stamina.x[t.index(max(t))],
                                        Stamina.y[t.index(max(t))]), 
                                5e3, color="orange", fill=True))
        plt.show()
    graph_EarthSurface()

    def graph_EarthMoon_Exact():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Earth Moon')
        plt.plot(Moon.x, Moon.y, color="#377EB8")
        plt.plot(Stamina.x, Stamina.y, color="orange")
        plt.xlim(xmin=-Moon.dE)
        plt.xlim(xmax=Moon.dE)
        plt.ylim(ymin=-Moon.dE)
        plt.ylim(ymax=Moon.dE)
        plt.hlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        plt.vlines(0, -Moon.dE, Moon.dE, "black", linewidth=1.0)
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), 
                                Moon.radius, color="#377EB8", fill=None))
        ax.add_patch(plt.Circle((Stamina.y[t.index(max(t))],
                                        Stamina.y[t.index(max(t))]), 
                                5e4, color="orange", fill=True))
        plt.show()
    graph_EarthMoon_Exact()

    def graph_EarthMoon_Margin():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Earth Moon')
        plt.plot(Moon.x, Moon.y, color="#377EB8")
        plt.plot(Stamina.x, Stamina.y, color="orange")
        plt.xlim(xmin=-Moon.dE - 1e8)
        plt.xlim(xmax=Moon.dE + 1e8)
        plt.ylim(ymin=-Moon.dE - 1e8)
        plt.ylim(ymax=Moon.dE + 1e8)
        plt.hlines(0, -Moon.dE - 1e8, Moon.dE + 1e8, "black", linewidth=1.0)
        plt.vlines(0, -Moon.dE - 1e8, Moon.dE + 1e8, "black", linewidth=1.0)
        plt.hlines(Moon.dE, -Moon.dE, Moon.dE, "gray", linewidth=0.5)
        plt.hlines(-Moon.dE, -Moon.dE, Moon.dE, "gray", linewidth=0.5)
        plt.vlines(Moon.dE, -Moon.dE, Moon.dE, "gray", linewidth=0.5)
        plt.vlines(-Moon.dE, -Moon.dE, Moon.dE, "gray", linewidth=0.5)
        ax.add_patch(plt.Circle((0,0), Moon.dE, color="gray", fill=None))
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), 
                                Moon.radius, color="#377EB8", fill=None))
        ax.add_patch(plt.Circle((Stamina.y[t.index(max(t))],
                                        Stamina.y[t.index(max(t))]), 
                                5e4, color="orange", fill=True))
        plt.show()
    graph_EarthMoon_Margin()

    def graph_CrashSite():
        fig=plt.figure(figsize=(8,8), dpi = 300)
        ax=fig.add_subplot(1,1,1)
        plt.title('Crash Site')
        plt.plot(Stamina.x, Stamina.y, color="orange")
        plt.plot(Stamina.x[i], Stamina.y[i], color="orange")
        plt.xlim(xmin=Stamina.x[i]-1.0e5)
        plt.xlim(xmax=Stamina.x[i]+1.0e5)
        plt.ylim(ymin=Stamina.y[i]-1.0e5)
        plt.ylim(ymax=Stamina.y[i]+1.0e5)
        ax.add_patch(plt.Circle((0,0), Earth.radius, color="green", fill=None))
        ax.add_patch(plt.Circle((Stamina.x[t.index(max(t))],Stamina.y[t.index(max(t))]), 1e3, color="orange", fill=True))
        ax.add_patch(plt.Circle((Moon.x[t.index(max(t))],Moon.y[t.index(max(t))]), Moon.radius, color="#377EB8", fill=None))
        plt.show()

    if Stamina.crash:
        graph_CrashSite()
        print("Stamina crashed at", t[i+1], "s")
    else:
        print("Stamina didn't crash.")

if __name__ == "__main__":
    main()
