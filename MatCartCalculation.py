import math, numpy as np



def ellipsoid(elips_name):
    if elips_name=='WGS84':
        a=6378137
        b=6356752.314
    elif elips_name=='Krasovskiy':
        a = 6378245
        b = 6356863.019
    else:
        a=1
        b=1
    e2=1-(b/a)**2
    ee2=(a/b)**2-1
    n=(a-b)/(a+b)
    m=(a*a-b*b)/(a*a+b*b)

    return a,b,e2

def NFirstVertical(B_rad,elips_name):
    #Радиус кривызны первого вертикала
    a, b, e2 = ellipsoid(elips_name)
    W=math.sqrt(1-e2*math.sin(B_rad)*math.sin(B_rad))
    N = a / W
    return  N
def rParalell(B_rad,elips_name):
    #Радиус кривизны параллели
    a, b, e2 = ellipsoid(elips_name)
    W = math.sqrt(1 - e2 * math.sin(B_rad) * math.sin(B_rad))
    r=a*math.cos(B_rad)/W
    return  r
def MMeridian(B_rad,elips_name):
    #Радиус кривызны меридиана
    a, b, e2 = ellipsoid(elips_name)
    W = math.sqrt(1 - e2 * math.sin(B_rad) * math.sin(B_rad))
    M = a * (1 - e2) / W**3
    return  M
def Geografic2CartesianXYZ(B_rad,L_rad,h_ell,elips_name):
    a, b, e2 = ellipsoid(elips_name)
    Nn=NFirstVertical(B_rad,elips_name)
    X=(Nn+h_ell)*math.cos(B_rad)*math.cos(L_rad)
    Y=(Nn+h_ell)*math.cos(B_rad)*math.sin(L_rad)
    Z=(Nn*(b/a)**2+h_ell)*math.sin(B_rad)
    return X,Y,Z
def X_DugaMeridiana2(B_rad,elips_name):
    a, b, e2 = ellipsoid(elips_name)
    A0=a*(1-e2)
    C0=A0*(1+e2*3/4+e2*e2*45/64+e2*e2*e2*175/256)
    C2=A0*(e2*3/8+e2*e2*15/32+e2*e2*e2*525/1024)
    C4=A0*(e2*e2*15/256+e2*e2*e2*105/1024)
    C6=A0*(e2*e2*e2*35/3072)
    print(C0, C2, C4,C6)
    X=C0*B_rad-C2*math.sin(2*B_rad)+C4*math.sin(4*B_rad)-C6*math.sin(6*B_rad)
    return X
def X_DugaMeridiana(B_rad,elips_name):
    a, b, e2 = ellipsoid(elips_name)
    s0=a*(1-e2)
    s2=e2*s0*3/2
    s4=e2*s2*5/4
    s6=e2*s4*7/6
    s8=e2*s6*9/8
    C0=s0+s2/2+s4*3/8+s6*5/16+s8*35/128
    C2=1/4*(s2+s4+s6*15/16+s8*7/8)
    C4=1/32*(s4+s6*3/2+s8*7/4)
    C6=1/96*(s6*1/2+s8)

    X=C0*B_rad-C2*math.sin(2*B_rad)+C4*math.sin(4*B_rad)-C6*math.sin(6*B_rad)
    return X
def B_Shirota_rad(X_DugaMeridiana,elips_name):
    a, b, e2 = ellipsoid(elips_name)
    s0=a*(1-e2)
    s2=e2*s0*3/2
    s4=e2*s2*5/4
    s6=e2*s4*7/6
    s8=e2*s6*9/8
    C0=s0+s2/2+s4*3/8+s6*5/16+s8*35/128
    C2=1/4*(s2+s4+s6*15/16+s8*7/8)
    C4=1/32*(s4+s6*3/2+s8*7/4)
    C6=1/96*(s6*1/2+s8)

    Betta=X_DugaMeridiana/C0
    D2=C2/C0*(1+C4/C0-C2*C2/(2*C0*C0))
    D4=C2*C2/(C0*C0)-C4/C0
    D6=C6/C0-3*C2/C0*(C4/C0-C2*C2/(2*C0*C0))
    print( D2*10**10, D4*10**10, D6*10**10)
    B=Betta+D2*math.sin(2*Betta)+D4*math.sin(4*Betta)+D6*math.sin(6*Betta)
    return B
def ENU2Cartesian(XENU, YENU, ZENU, Fi0_rad, Lym0_rad, Height0, elips_name):
    X0,Y0,Z0=Geografic2CartesianXYZ(Fi0_rad, Lym0_rad, Height0, elips_name)
    ENU_XYZ = np.array([[XENU], [YENU], [ZENU]], float)  # вертикальный массив
    ENURotation = np.zeros((3, 3), float)
    ########MakeRotationMatrix#######
    cos_fi = math.cos(Fi0_rad)
    cos_lym = math.cos(Lym0_rad)
    sin_fi = math.sin(Fi0_rad)
    sin_lym = math.sin(Lym0_rad)
    #
    ENURotation[0, 0] = - sin_lym
    ENURotation[0, 1] = - sin_fi * cos_lym
    ENURotation[0, 2] = cos_fi * cos_lym
    #
    ENURotation[1, 0] = cos_lym
    ENURotation[1, 1] = - sin_fi * sin_lym
    ENURotation[1, 2] = cos_fi * sin_lym
    #
    ENURotation[2, 0] = 0
    ENURotation[2, 1] = cos_fi
    ENURotation[2, 2] = sin_fi
    #print('ENURotation\n', ENURotation)
    deltaXCartesian, deltaYCartesian, deltaZCartesian = np.dot(ENURotation, ENU_XYZ)  # М
    LidarXCartesian =X0+ deltaXCartesian
    LidarYCartesian = Y0 + deltaYCartesian
    LidarZCartesian = Z0 + deltaZCartesian
    return LidarXCartesian, LidarYCartesian, LidarZCartesian

def Cartesian2Geografic(XCartesian, YCartesian, ZCartesian, elips_name, delta=0.5 * 10 ** -10):
    def fi_rad(XCartesian, YCartesian,ZCartesian,S):
        XY = (XCartesian * XCartesian + YCartesian * YCartesian)
        if XY !=0:
            fi=math.asin(S)
        else:
            fi=math.pi/2*ZCartesian/abs(ZCartesian)
        return fi
    def lym_rad(LidarXCartesian, LidarYCartesian):
        XYdist=math.sqrt(LidarXCartesian*LidarXCartesian + LidarYCartesian*LidarYCartesian)
        if LidarXCartesian>=0:
            lym=math.asin(LidarYCartesian/XYdist)
        elif LidarXCartesian < 0 and LidarYCartesian >=0:
            lym = math.acos(LidarXCartesian / XYdist)
        elif LidarXCartesian < 0 and LidarYCartesian <0:
            lym = -math.acos(LidarXCartesian / XYdist)
        else:
            lym=0
        return lym
    diff=1
    S1=0
    a, b, e2 = ellipsoid(elips_name)
    while diff>delta:
        N=a/math.sqrt((1-e2*S1*S1))
        P=e2*N*S1
        Q=math.sqrt(XCartesian * XCartesian + YCartesian * YCartesian + (ZCartesian + P) * (ZCartesian + P))
        S2= (ZCartesian + P) / Q
        diff=abs(S2-S1)
        S1=S2
        #print(diff)
    Fi_rad=fi_rad(XCartesian, YCartesian, ZCartesian, S2)
    Lym_rad=lym_rad(XCartesian, YCartesian)
    Height = Q-N
    return Fi_rad, Lym_rad, Height

def XYGaussKruger(B_rad, L_rad, L0_rad, m0=0.996,elips_name='WGS84' ):
    a, b, e2 = ellipsoid(elips_name)
    ee2 = (a / b) ** 2 - 1
    sin_B=math.sin(B_rad)
    cos_B=math.cos(B_rad)
    t=math.tan(B_rad)
    nu2=ee2*cos_B*cos_B
    dl=L_rad-L0_rad
    N=NFirstVertical(B_rad,elips_name)
    X=X_DugaMeridiana(B_rad,elips_name)
    ###
    XX1=5 - t*t + 9*nu2 + 4*nu2*nu2
    XX2=61-58*t*t + t**4 + 270*nu2 - 330*nu2*t*t
    x=m0*(X+N*sin_B*(1/2*cos_B*dl*dl + 1/24*cos_B**3*XX1*dl**4 + 1/720*cos_B**5*XX2*dl**6))
    YY1=1 - t*t + nu2
    YY2=5 - 18*t*t + t**4 + 14*nu2 - 58*nu2*t*t
    y=m0*N*cos_B*(dl+1/6*cos_B**2*YY1*dl**3 + 1/120*cos_B**4*YY2*dl**5)
    return x,y

def find_Cenrtal_Meridian(lym1, fi1,lym2,fi2,Y1,X1,Y2,X2):
    print('start')
    for L0_min in range(180*60):
        #L0_rad = 66
        L0_rad = math.radians(L0_min/60)
        
        X1calc,Y1calc=XYGaussKruger(fi1, lym1, L0_rad, 1,'Krasovskiy' )
        #print(X1calc,Y1calc)
        X2calc,Y2calc=XYGaussKruger(fi2, lym2, L0_rad, 1,'Krasovskiy' )
        dX1=X1calc-X1
        dY1=Y1calc-Y1
        dX2=X2calc-X2
        dY2=Y2calc-Y2
        dX=abs(dX1-dX2)
        dY=abs(dY1-dY2)
        delta=math.sqrt(dX*dX+dY*dY)
        if L0_min%60==0:
            pass
            #print(math.degrees(L0_rad),dX1,dX2)
            #print(math.degrees(L0_rad),dX,dY)
        if delta<30:
            print(delta,math.degrees(L0_rad),(dX1+dX2)/2,(dY1+dY2)/2)

#############################################

if __name__=="__main__":
    ################################LIDAR SAMPLE DATA##########################
    LidarAzimut_rad = math.radians(144.875)
    LidarVertical_rad = math.radians(-0.667)
    LidarDistance = 163.852
    #########
    YawLidar_rad = math.radians(0)
    PitchLidar_rad = math.radians(90)
    RollLidar_rad = math.radians(0)
    # Смещения измерядтся по направлению осей координат самолета от IMU  и прибавляются к коордиантам Лидара
    XshiftLidar = 0.1
    YshiftLidar = 0
    ZshiftLidar = 0.05
    ######
    B_rad = math.radians(45.10045587)
    L_rad = math.radians(38.24548502)
    h_ell = 166.1319532
    YawVehicle_rad = math.radians(254.607332)
    PitchVehicle_rad = math.radians(3)
    RollVehicle_rad = math.radians(0.07332)
    elips_name = 'WGS84'

    ###########################################################################
    #print('B', math.degrees(B_rad))
    #X=X_DugaMeridiana(B_rad,elips_name)
    #print('S_DugaMeridiana(B_rad,elips_name)', X)
    #B=B_Shirota_rad(X,elips_name)
    #print('B_Shirota_rad(X_DugaMeridiana,elips_name)',math.degrees(B) )
    #print('XYGaussKruger(B_rad, L_rad, L0_rad, m0=1 )',XYGaussKruger(B_rad, L_rad, math.radians(39)) )
    #XY_etalon = (84137.1960142879, 4541135.27680016)
    #print('XYGaussKruger(B_rad, L_rad, L0_rad, m0=1 )',XYGaussKruger(math.radians(41), math.radians(38), math.radians(37), 1,'Krasovskiy' ) )
    # 67.376959733,70.420666317,67.48063731619,70.40724066520
    # 251658.1338,309932.6426,255489.2306,308317.1073
    # 67.38450,70.55152,67.18689,70.86150
    # 252310.8462,324505.0575,245997.440,359284.413
    find_Cenrtal_Meridian(math.radians(67.38450),math.radians(70.55152),math.radians(67.18689),math.radians(70.86150),252310.8462,324505.0575,245997.440,359284.413)
    #find_Cenrtal_Meridian(math.radians(67.41119859),math.radians(70.66029005),math.radians(67.4401892),math.radians(70.4069995),253692.13465,336684.05248,254146.95552825,308382.32057980)
