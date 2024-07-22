import freecad
import Part
# from freecad import Base
import numpy as np
import math
from itertools import product
import pathlib as Path
cwd = Path.Path().resolve()
path = str(cwd.parents[0])
import sys
sys.path.append(path) 

# CONSTANT PARAMETERS
APERTURE_AREA = 8000 * 8000 # mm2
TOWER_GAP = 1800 # mm
NUMBER_MIRRORS = 4 ** 2
SPHERE_DIAMETER = 300 # mm
MIRRORS_GAP = 20 # mm
Primary_Glass_thickness = 4 # mm
Secondary_Glass_thickness = 4 # mm
CPC_Glass_thickness = 4 # mm
Frame_Width = 70 # mm

# VARIABLE PARAMETERS
P_FOCAL_DISTANCE = np.arange(1000, 4501, 500)
HYPERPOLA_FOCUS = np.arange(100, 901,200)
P_SPHERE = np.arange(-400, 601, 200)
APERTURE_SPHERE = np.arange(200, 241, 20)
ACCEPTANCE_ANGLE = np.arange(14, 24.5, 2)
TRUNCATION_FACTOR = np.arange(0.20, 0.81, 0.2)

parameters_set = product(P_FOCAL_DISTANCE, HYPERPOLA_FOCUS, P_SPHERE,
                        APERTURE_SPHERE, ACCEPTANCE_ANGLE, TRUNCATION_FACTOR)
    
def single_design(arg):
    (
        P_FOCAL_DISTANCE,
        HYPERPOLA_FOCUS,
        P_SPHERE,
        APERTURE_SPHERE,
        ACCEPTANCE_ANGLE,
        TRUNCATION_FACTOR,
    ) = arg
    # Alias for active Document
    doc = App.newDocument()
    
    parabola_size = ( 2 * ( APERTURE_AREA ** 0.5 + TOWER_GAP + MIRRORS_GAP * (NUMBER_MIRRORS ** 0.5 - 1) ) ** 2 ) ** 0.5
    h_parabola = (parabola_size/2)**2 / (4 * P_FOCAL_DISTANCE)
    # create a Part object that is a Parabola in the XY plane (the parabola is infinite).
    parabola_curve = Part.Parabola()
    # define de Focal distance in the X axe
    parabola_curve.Focal = P_FOCAL_DISTANCE
    # create an edge from the parabola curve with the limits in the Y axe
    parabola_size = ( 2 * ( APERTURE_AREA ** 0.5 + TOWER_GAP + MIRRORS_GAP * (parabola_size ** 0.5 - 1) ) ** 2 ) ** 0.5
    edge_parabola = Part.Edge(parabola_curve, parabola_size / 2, 0)
    # adds a Part object type to the document and assigns the shape representation of the edge_p
    Parabola = doc.addObject("Part::Feature","Parabola")
    Parabola.Shape = edge_parabola
    # making a closed wire in order to revolute
    p1 = Part.Vertex(0, 0, 0)
    p2 = Part.Vertex(edge_parabola.BoundBox.XMax, 0, 0)
    p3 = Part.Vertex(edge_parabola.BoundBox.XMax, edge_parabola.BoundBox.YMax, 0)
    Line1 = Part.Edge(p1, p2)
    Line2 = Part.Edge(p2, p3)
    wire_parabola = Part.Wire([Line1, Line2, edge_parabola])
    face_parabola = Part.Face(wire_parabola)
    face_parabola.rotate(FreeCAD.Vector(0.0,0.0,0.0),FreeCAD.Vector(0.0,1.0,0.0),-90.0)
    face_parab = doc.addObject("Part::Feature", "face_parab")
    face_parab.Shape = face_parabola
    doc.recompute()

    # center mirrors
    center_mirrors = {}
    X_mirrors = {}
    Y_mirrors = {}
    projected_mirror_length = (APERTURE_AREA / NUMBER_MIRRORS) ** 0.5
    X_mirrors[1] = TOWER_GAP / 2 + projected_mirror_length / 2
    Y_mirrors[1] = projected_mirror_length / 2 + MIRRORS_GAP / 2
    
    for i in range(2, int(NUMBER_MIRRORS ** 0.25) + 1, 1):
    	X_mirrors[i] = X_mirrors[i - 1] + MIRRORS_GAP + projected_mirror_length
    	Y_mirrors[i] = Y_mirrors[i - 1] + MIRRORS_GAP + projected_mirror_length
    
    X_mirrors[-1] = - TOWER_GAP / 2 - projected_mirror_length / 2
    Y_mirrors[-1] = - projected_mirror_length / 2 - MIRRORS_GAP / 2
    
    for i in range(-2, -int(NUMBER_MIRRORS ** 0.25) - 1, -1):
    	X_mirrors[i] = X_mirrors[i + 1] - MIRRORS_GAP - projected_mirror_length
    	Y_mirrors[i] = Y_mirrors[i + 1] - MIRRORS_GAP - projected_mirror_length
    
    k = 1
    
    for i in range(-int(NUMBER_MIRRORS ** 0.25), int(NUMBER_MIRRORS ** 0.25) + 1):
    	for j in range(-int(NUMBER_MIRRORS ** 0.25), int(NUMBER_MIRRORS ** 0.25) + 1):
    		if i != 0 and j != 0:
    			center_mirrors[k] = App.Base.Vector(X_mirrors[i], Y_mirrors[j], 0)
    			k = k + 1
    center_mirrors[k] = App.Base.Vector(0, - Y_mirrors[j], 0)
    
    
    for k in range(1, NUMBER_MIRRORS + 2, 1):
        center_mirror = center_mirrors[k]
        box_width = projected_mirror_length
        box_Length = projected_mirror_length
        if k > NUMBER_MIRRORS:
            box_Length = TOWER_GAP -  2 * MIRRORS_GAP 
        box_Height = doc.face_parab.Shape.BoundBox.ZMax + 10
        actual_box = doc.addObject("Part::Box","actual_box")
        actual_box.Width = box_width
        actual_box.Length = box_Length
        actual_box.Height = box_Height
        center_vertex = App.Base.Vector(center_mirror[0] - box_Length / 2, center_mirror[1] - box_width / 2, -10)
        actual_box.Placement= App.Base.Placement(center_vertex,App.Rotation(App.Vector(0,0,1),0))
        # making a revolution of the face parabola in order to contruct de paraboloid
        # making a revolution of the parabola edge to construct the paraboloid
        Label_paraboloid_i = 'Paraboloid_{0}'.format(int(abs(k)))
        doc.addObject("Part::Revolution", Label_paraboloid_i)
        doc.ActiveObject.Source = doc.face_parab
        doc.ActiveObject.Axis = (0.0, 0.0, 1.0)
        doc.ActiveObject.Base = (0.0, 0.0, 0.0)
        doc.ActiveObject.Angle = 360.0
        doc.ActiveObject.Solid = True
        doc.ActiveObject.AxisLink = None
        doc.ActiveObject.Symmetric = False
        doc.recompute()
        # creating the mirror shape with common tool
        doc.addObject("Part::MultiCommon","common")
        doc.common.Shapes = [doc.actual_box,doc.getObject(Label_paraboloid_i)]
        doc.recompute()
        # this is the segment paraboloid face for the heliostat well positioned
        face = doc.common.Shape.Faces[4]
        if k == 16:
            face16 = face
        Label_mirror_i = 'Mirror_{0}'.format(int(abs(k)))
        actual_mirror = doc.addObject("Part::Feature","Label_mirror_i")
        doc.ActiveObject.Label = Label_mirror_i.__add__("(HRC_1)")
        doc.ActiveObject.Shape = face
        doc.removeObject(Label_paraboloid_i)
        doc.addObject('PartDesign::Body','Body')
        doc.addObject('PartDesign::FeatureBase','Clone')
        doc.getObject('Clone').BaseFeature = actual_mirror
        doc.getObject('Clone').Placement = actual_mirror.Placement
        doc.getObject('Clone').setEditorMode('Placement',0)
        doc.getObject('Body').Group = [doc.getObject('Clone')]
        doc.getObject('Body').Tip = doc.getObject('Clone')
        doc.getObject('Body').newObject('PartDesign::Pad','Pad')
        doc.getObject('Pad').Profile = doc.getObject('Clone')
        doc.getObject('Pad').Length = Primary_Glass_thickness
        doc.getObject('Pad').Reversed = 1
        doc.recompute()
        doc.recompute()
        Label_glass_i = 'Glass_{0}'.format(int(abs(k)))
        doc.addObject('Part::Feature',Label_glass_i)
        doc.ActiveObject.Shape = doc.Pad.Shape
        # doc.ActiveObject.Label= Label_glass_i.__add__("(SiO2)")
        doc.getObject('Body').removeObjectsFromDocument()
        doc.removeObject('common')
        doc.removeObject('actual_box')
        doc.removeObject('Body')
    
    doc.removeObject('Parabola')
    doc.removeObject('face_parab')
    
    
    XX = face16.BoundBox.YMax
    YY = face16.BoundBox.XMax
    ZZ = face16.BoundBox.ZMax
    D_half = (XX**2 + YY**2 ) **0.5
    _angle = np.arctan(P_FOCAL_DISTANCE/D_half)
    # ray deviated due to optical errors
    s_sun = 4.65/2/1000
    s_m = 2.4529/2/1000
    s_s = 2/1000
    s_t = 0.2 * np.pi / 180
    s_a = 2/1000
    sigma = ((2*s_sun)**2 + (2*s_m)**2 + (2*s_s)**2 + (2*s_t)**2 + (2*s_a)**2)**0.5
    _angle = _angle + sigma
    Z_Point = np.tan(_angle) * D_half
    WIDTH_SECONDARY = 2 * (Z_Point - (P_FOCAL_DISTANCE - HYPERPOLA_FOCUS)) * D_half / Z_Point / 2 ** 0.5 ## 
    # WIDTH_SECONDARY = TOWER_GAP - 2 * Frame_Width

    ## CPC receiver
    # Inputs for the CPC
    a_width = APERTURE_SPHERE # absorber width
    angle_ = np.arcsin((a_width/2) / (SPHERE_DIAMETER/2))
    height_absorber = P_SPHERE +  (SPHERE_DIAMETER / 2) *np.cos(angle_) #  thermal absorber height from the primary mirrors
    acceptance_angle = ACCEPTANCE_ANGLE # acceptance angle for the CPC
    theta_ = acceptance_angle * np.pi / 180 # acceptance angle for the CPC in radians
    
    # Parameters for geometry determination
    a = a_width / 2 # half absorber width
    aPrima = a / math.sin(theta_) # half absorber width
    focal_distance_CPC = a * ( 1 + math.sin(theta_)) # CPC focal distance
    h_CPC = focal_distance_CPC * math.cos(theta_) / (math.sin(theta_)**2) # CPC height without truncation  
    h_CPC_T = h_CPC * TRUNCATION_FACTOR #2 # CPC height truncated
    H_CPC_aperture = height_absorber - h_CPC_T # height of the CPC aperture from primary mirrors
    b = a / math.tan(theta_) # Intersection between acceptance lines
    
    # Algortihm to find the phi angle coordinates according to the CPC truncation.
    # The phi angle is the aperture angle of the truncated CPC.
    for phi_ in np.arange(theta_, np.pi, 1E-6):
    	hi = focal_distance_CPC * math.cos(phi_ - theta_) / (math.sin(phi_ / 2)**2)
    	if (hi < h_CPC_T):
    		break
    
    a_CPC_T = focal_distance_CPC * math.sin(phi_ - theta_)/((math.sin(phi_ / 2)**2)) - a # half CPC truncated aperture 
    # m and n parameters for the CPC parabolas
    n = 2 * a * math.cos(theta_) 
    m = (a + a_CPC_T) * math.sin(phi_)/math.sin(phi_ - theta_)
    
    #--------------------------------------- CPC ------------------------------------------------------
    x_pos = a - focal_distance_CPC * math.sin(theta_) # X placement parabola for CPC
    y_pos = focal_distance_CPC * math.cos(theta_) # Y placement parabola for CPC
    #Create parabola CPC
    parab_CPC = Part.Parabola()
    # define de Focal distance in the X axe
    parab_CPC.Focal = focal_distance_CPC
    # create an edge from the parabola curve with the limits in the Y axe
    half_parab_CPC = Part.Edge(parab_CPC, m, n)
    
    # transformating the face_parabola and creating a face Part
    half_parab_CPC.rotate(App.Vector(0.0,0.0,0.0),App.Vector(0.0,0.0,1.0),-theta_ * 180/np.pi)
    half_parab_CPC.rotate(App.Vector(0,0,0),App.Vector(0,0,1),90)
    half_parab_CPC.rotate(App.Vector(0,0,0),App.Vector(1,0,0),90)
    
    half_parab_CPC.translate(App.Base.Vector(x_pos,0,-y_pos))
    half_parab_CPC.translate(App.Base.Vector(0,0,height_absorber))
    # Part.show(half_parab_CPC)
    # creating a face with half parabola form
    p1 = Part.Vertex(half_parab_CPC.BoundBox.XMax,0,half_parab_CPC.BoundBox.ZMin)
    p2 = Part.Vertex(0,0,half_parab_CPC.BoundBox.ZMin)
    p3 = Part.Vertex(0,0,half_parab_CPC.BoundBox.ZMax)
    p4 = Part.Vertex(half_parab_CPC.BoundBox.XMin,0,half_parab_CPC.BoundBox.ZMax)
    Line1 = Part.Edge(p1,p2)
    Line2 = Part.Edge(p2,p3)
    Line3 = Part.Edge(p3,p4)
    wire_parabola = Part.Wire([half_parab_CPC,Line1,Line2,Line3])
    face_parabola = Part.Face(wire_parabola) 
    face_2parab = doc.addObject("Part::Feature","face_2parab")
    face_2parab.Shape = face_parabola
    doc.recompute()
    # making a revolutionof the face Part
    doc.addObject("Part::Revolution","CPC_revolution")
    doc.CPC_revolution.Source = doc.face_2parab
    doc.CPC_revolution.Axis = (0.0,0.0,1.0)
    doc.CPC_revolution.Base = (0.0,0.0,0.0)
    doc.CPC_revolution.Angle = 360.0
    doc.CPC_revolution.Solid = True
    doc.CPC_revolution.AxisLink = None
    doc.CPC_revolution.Symmetric = False
    # doc.face_2parab.Visibility = False
    doc.recompute()
    
    # this is the paraboloid face for the heliostat
    face_CPC = doc.CPC_revolution.Shape.Faces[0]
    
    # creating a Part for the specular surface Mirror
    mirror_CPC = doc.addObject("Part::Feature",'mirror_CPC')
    doc.mirror_CPC.Label='mirror_CPC(HRC_1)'
    doc.ActiveObject.Shape=face_CPC
    # doc.ActiveObject.Placement=Base.Placement(Base.Vector(0,0,0),freecad.Rotation(freecad.Vector(0,0,1),0))
    # gui.ActiveObject.ShapeColor = (0.0,1.0,1.0)
    
    doc.addObject('PartDesign::Body','Body')
    doc.addObject('PartDesign::FeatureBase','Clone')
    doc.getObject('Clone').BaseFeature = doc.getObject('mirror_CPC')
    doc.getObject('Clone').Placement = doc.getObject('mirror_CPC').Placement
    doc.getObject('Clone').setEditorMode('Placement',0)
    doc.getObject('Body').Group = [doc.getObject('Clone')]
    doc.getObject('Body').Tip = doc.getObject('Clone')
    
    doc.getObject('Body').newObject('PartDesign::Pad','Pad')
    doc.getObject('Pad').Profile = doc.getObject('Clone')
    doc.getObject('Pad').Length = CPC_Glass_thickness
    doc.recompute()
    doc.recompute()
    
    doc.addObject('Part::Feature','CPC_Glass')
    doc.CPC_Glass.Shape = doc.Pad.Shape
   
    # doc.CPC_Glass.Label='CPC_Glass(SiO2)'
    
    doc.getObject('Body').removeObjectsFromDocument()
    doc.removeObject('Body')
       
    ## Sphere for the receiver
    shift = ((SPHERE_DIAMETER/2)**2 - (a_width/2)**2)**0.5
    #create a FreeCAD object with Cylinder attributes
    sphere = doc.addObject("Part::Sphere","sphere")
    sphere.Radius = SPHERE_DIAMETER/2.0
    sphere.Angle3 = 360.0
    sphere.Placement = App.Base.Placement(App.Vector(0.0,0.0,P_SPHERE),App.Rotation(App.Vector(1,0.0,0.0),0))
    doc.recompute()
        
    # creating the sphere shape with Cut tool
    doc.addObject("Part::Cut","Cut")
    doc.Cut.Base = doc.sphere
    doc.Cut.Tool = doc.CPC_revolution
    # gui.sphere.Visibility=False
    # gui.CPC_revolution.Visibility=False
    
    doc.recompute()
    # this is the paraboloid face for the heliostat
    face_abs = doc.Cut.Shape.Faces[0]
    # creating a Part for the Receiver
    doc.addObject("Part::Feature",'Receiver')
    doc.Receiver.Label='Receiver(MoSi2Si3N4)'
    doc.Receiver.Shape=face_abs
    # doc.ActiveObject.Placement=Base.Placement(Base.Vector(0,0,0),freecad.Rotation(freecad.Vector(0,0,1),0))
    # gui.ActiveObject.ShapeColor = (0.0,1.0,1.0)
    
    doc.removeObject('CPC_revolution')
    doc.removeObject('sphere')
    doc.removeObject('Cut')
    doc.removeObject('face_2parab')
    
    ## Secondary Reflector
    main_focus_distance = P_FOCAL_DISTANCE - height_absorber
    two_c = main_focus_distance  # this is the focal distance
    c = two_c / 2 
    two_a = two_c - 2 * HYPERPOLA_FOCUS # this is two times the major radius of the hyperbola
    a = two_a / 2
    print(a, P_FOCAL_DISTANCE, height_absorber, HYPERPOLA_FOCUS)
    if a < 0:
       return None
    b = (c**2 - a **2)**0.5 # this is the minor radius of the parabola
    center = App.Vector(height_absorber + c,0,0) # this is the vertex of the hyperbola
    
    endpoint_1 = - (np.arcsinh(WIDTH_SECONDARY / b))
    endpoint_2 = (np.arcsinh(WIDTH_SECONDARY / b))
    endpoint_1 = 0
    
    # create a Part object that is a Hyperbola in the XY plane (the hyperbola is infinite).
    hyperbola_curve = Part.Hyperbola()
    # define the Major Radius in the X axe
    hyperbola_curve.MajorRadius = a
    hyperbola_curve.MinorRadius = b
    hyperbola_curve.Center = center
    # create an edge from the parabola curve with the limits in the Y axe
    arc_hyperbola = Part.ArcOfHyperbola(hyperbola_curve, endpoint_1, endpoint_2)
    edge_hyperbola = Part.Edge(arc_hyperbola)
    #Part.show(edge_hyperbola)
    p1 = Part.Vertex(edge_hyperbola.BoundBox.XMax,edge_hyperbola.BoundBox.YMax,0)
    p2 = Part.Vertex(edge_hyperbola.BoundBox.XMax,0,0)
    p3 = Part.Vertex(edge_hyperbola.BoundBox.XMin,edge_hyperbola.BoundBox.YMin,0)
    Line1 = Part.Edge(p1,p2)
    Line2 = Part.Edge(p2,p3)
    wire_hyperbola = Part.Wire([edge_hyperbola,Line1,Line2])
    face_hyperbola = Part.Face(wire_hyperbola)
    
    # transformating the face_hyperbola and creating a face Part
    face_hyperbola.rotate(App.Vector(0.0,0.0,0.0),App.Vector(0.0,1.0,0.0),-90.0)
    face_hyp = doc.addObject("Part::Feature","face_hyp")
    face_hyp.Shape = face_hyperbola
    # making a revolutionof the face Part
    doc.addObject("Part::Revolution","revolve_hyperbola")
    doc.revolve_hyperbola.Source = doc.face_hyp
    doc.revolve_hyperbola.Axis = (0.0,0.0,1.0)
    doc.revolve_hyperbola.Base = (0.0,0.0,0.0)
    doc.revolve_hyperbola.Angle = 360.0
    doc.revolve_hyperbola.Solid = True
    doc.revolve_hyperbola.AxisLink = None
    doc.revolve_hyperbola.Symmetric = False
    
    doc.recompute
        
    Width = WIDTH_SECONDARY
    Length = WIDTH_SECONDARY
    Height = WIDTH_SECONDARY * 2 
    box = doc.addObject("Part::Box","box")
    box.Width = Width
    box.Length = Length
    box.Height = Height
    box.Placement=App.Base.Placement(App.Base.Vector(-Length/2,-Width / 2,height_absorber + two_c - HYPERPOLA_FOCUS - 10),App.Base.Rotation(App.Base.Vector(0,0,1),0))
    
    doc.addObject("Part::MultiCommon","common")
    doc.common.Shapes = [doc.box,doc.revolve_hyperbola]
    
    doc.recompute()
    # this is the paraboloid face for the secondary
    face = doc.common.Shape.Faces[4]
    #creating a Part of the paraboloid face
    doc.addObject('Part::Feature','Secondary')
    doc.Secondary.Shape = face
    doc.Secondary.Label='Secondary(HRC_1)'
    Z_elements = doc.Secondary.Shape.BoundBox.ZMax
    
    doc.addObject('PartDesign::Body','Body')
    doc.addObject('PartDesign::FeatureBase','Clone')
    doc.getObject('Clone').BaseFeature = doc.getObject('Secondary')
    doc.getObject('Clone').Placement = doc.getObject('Secondary').Placement
    doc.getObject('Clone').setEditorMode('Placement',0)
    doc.getObject('Body').Group = [doc.getObject('Clone')]
    doc.getObject('Body').Tip = doc.getObject('Clone')
    
    doc.getObject('Body').newObject('PartDesign::Pad','Pad')
    doc.getObject('Pad').Profile = doc.getObject('Clone')
    doc.getObject('Pad').Length = Secondary_Glass_thickness
    doc.recompute()
    doc.recompute()
    
    doc.addObject('Part::Feature','Secondary_Glass')
    doc.Secondary_Glass.Shape = doc.Pad.Shape
    # doc.Secondary_Glass.Label='Secondary_Glass(SiO2)'
    
    doc.getObject('Body').removeObjectsFromDocument()
    doc.removeObject('Body')
    
    # to recalculate the whole document
    doc.recompute()
    
    doc.removeObject('face_hyp')
    doc.removeObject('box')
    doc.removeObject('common')
    doc.removeObject('revolve_hyperbola')
    doc.recompute()
    
    tower_width = TOWER_GAP - 2 * Frame_Width
    Width = tower_width
    Length = tower_width
    Height = 5000 
    box_1 = doc.addObject("Part::Box","box_1")
    box_1.Width = Width
    box_1.Length = Length
    box_1.Height = Height
    box_1.Placement=App.Base.Placement(App.Base.Vector(-Length/2,-Width / 2, P_SPHERE + SPHERE_DIAMETER - Height),App.Base.Rotation(App.Base.Vector(0,0,1),0))
    
    Width = tower_width
    Length = tower_width - SPHERE_DIAMETER * 1.5
    Height = SPHERE_DIAMETER * 2 
    box_2 = doc.addObject("Part::Box","box_2")
    box_2.Width = Width
    box_2.Length = Length
    box_2.Height = Height
    box_2.Placement=App.Base.Placement(App.Base.Vector(-Length/2,-Width / 2, P_SPHERE + SPHERE_DIAMETER - Height),App.Base.Rotation(App.Base.Vector(0,0,1),0))
    
    doc.addObject("Part::Cut","Tower")
    doc.Tower.Base = doc.box_1
    doc.Tower.Tool = doc.box_2
    doc.Tower.Label = "Tower(Opaque)"
    doc.recompute()
    
    Radius = SPHERE_DIAMETER / 2 * 1.3
    Height = tower_width - SPHERE_DIAMETER * 1.5
    cylinder = doc.addObject("Part::Cylinder","cylinder")
    cylinder.Radius = Radius
    cylinder.Height = Height
    # exchanger.Label = "Exchanger(Opaque)"
    cylinder.Placement=App.Base.Placement(App.Base.Vector(Height/2,0, P_SPHERE),App.Base.Rotation(App.Base.Vector(0,1,0),-90))
    doc.recompute()

    doc.addObject("Part::Cut","Cut")
    doc.Cut.Base = doc.cylinder
    doc.Cut.Tool = doc.CPC_Glass
    doc.recompute()

    # this is the exchanger face
    face_exchanger = doc.Cut.Shape.Faces[1]
    # creating a Part for the Receiver
    doc.addObject("Part::Feature",'Exchanger')
    doc.Exchanger.Label = 'Exchanger(Opaque)'
    doc.Exchanger.Shape = face_exchanger
    doc.recompute()
    doc.removeObject('cylinder')
    doc.removeObject('Cut')
    
    Width = Frame_Width
    Length = Frame_Width
    Height = Z_elements
    mast_1 = doc.addObject("Part::Box","Mast_1")
    mast_1.Width = Width
    mast_1.Length = Length
    mast_1.Height = Height
    mast_1.Label = "Mast_1(Opaque)"
    mast_1.Placement=App.Base.Placement(App.Base.Vector(tower_width / 2,-Length / 2, 0),App.Base.Rotation(App.Base.Vector(0,0,1),0))
    
    Width = Frame_Width
    Length = Frame_Width
    Height = Z_elements
    mast_2 = doc.addObject("Part::Box","Mast_2")
    mast_2.Width = Width
    mast_2.Length = Length
    mast_2.Height = Height
    mast_2.Label = "Mast_2(Opaque)"
    mast_2.Placement=App.Base.Placement(App.Base.Vector(- tower_width / 2  - Width,-Length / 2, 0),App.Base.Rotation(App.Base.Vector(0,0,1),0))
    
    Radius = Width / 2
    Height = tower_width
    mast_3 = doc.addObject("Part::Cylinder","Mast_3")
    mast_3.Radius = Radius
    mast_3.Height = Height
    mast_3.Label = "Mast_3(Opaque)"
    mast_3.Placement=App.Base.Placement(App.Base.Vector(tower_width / 2,0, Z_elements - 60 / 2),App.Base.Rotation(App.Base.Vector(0,1,0),-90))
    doc.recompute()
    a = list(arg)
    b = round(a[5], 2)
    a[5] = b
    print(a)
    if face_CPC.BoundBox.ZMax + 100 < face.BoundBox.ZMin:
        Label_drawing = "designs_set_1/design_{0}".format(a)+".FCStd"
        doc.saveAs(Label_drawing)

case = 1
# total_cases = len(list(parameters_set))

for arg in parameters_set:
    case = case + 1
    print(case / 17280 * 100, "%")
    actual_design = single_design(arg)
    
