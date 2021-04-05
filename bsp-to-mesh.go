package main

import (
    "github.com/galaco/bsp"
    "github.com/galaco/bsp/lumps"
    "github.com/galaco/bsp/primitives/face"
    "github.com/galaco/bsp/primitives/plane"
    "github.com/galaco/bsp/primitives/texinfo"
    "github.com/galaco/bsp/primitives/brush"
    "github.com/galaco/bsp/primitives/brushside"
    "github.com/galaco/bsp/primitives/model"
    "github.com/go-gl/mathgl/mgl32"
    "log"
    "math"
    "os"
    "encoding/json"
    "fmt"
    "strconv"
    "errors"
    "gonum.org/v1/gonum/stat/combin"
    "github.com/tchayen/triangolatte"
    "sort"
    "io"
    "github.com/golang-source-engine/stringtable"
    "encoding/binary"
    "strings"
    "bufio"
)

var bspId = 0
const brushCollisionsOnly = true

type BSPData struct {
    faces      []face.Face
    planes     []plane.Plane
    vertexes   []mgl32.Vec3
    surfEdges  []int32
    edges      [][2]uint16
    texInfos   []texinfo.TexInfo
    brushes    []brush.Brush
    brushSides []brushside.BrushSide
    models     []model.Model
}

func main() {
    if len(os.Args) < 2 {
        log.Fatal(fmt.Sprintf("Usage: %s <filename>", os.Args[0]))
    }

    filename := os.Args[1]
    file, err := bsp.ReadFromFile(filename)
    if err != nil {
        log.Fatal(err)
    }
    
    bspData := BSPData {
        faces:      file.Lump(bsp.LumpFaces).(*lumps.Face).GetData(),
        planes:     file.Lump(bsp.LumpPlanes).(*lumps.Planes).GetData(),
        vertexes:   file.Lump(bsp.LumpVertexes).(*lumps.Vertex).GetData(),
        surfEdges:  file.Lump(bsp.LumpSurfEdges).(*lumps.Surfedge).GetData(),
        edges:      file.Lump(bsp.LumpEdges).(*lumps.Edge).GetData(),
        texInfos:   file.Lump(bsp.LumpTexInfo).(*lumps.TexInfo).GetData(),
        brushes:    file.Lump(bsp.LumpBrushes).(*lumps.Brush).GetData(),
        brushSides: file.Lump(bsp.LumpBrushSides).(*lumps.BrushSide).GetData(),
        models:     file.Lump(bsp.LumpModels).(*lumps.Model).GetData(),
    }

    // load string table
    stringData := file.Lump(bsp.LumpTexDataStringData).(*lumps.TexDataStringData).GetData()
    stringTableArrayBytes, _ := file.Lump(bsp.LumpTexDataStringTable).Marshall()
    stringTableArray := make([]int32, len(stringTableArrayBytes)/4)
    for i := range stringTableArray {
        intBytes := stringTableArrayBytes[i*4:(i+1)*4]
        stringTableArray[i] = int32(binary.LittleEndian.Uint32(intBytes))
    }
    stringTable := stringtable.NewFromExistingStringTableData(stringData, stringTableArray)

    // entities
    entityLumpStringRaw := strings.Trim(file.Lump(bsp.LumpEntities).(*lumps.EntData).GetData(), "\x00")
    entityLumpStrings := strings.SplitAfter(entityLumpStringRaw, "}")
    entities := make([]map[string]string, 0)
    for _, entityLumpString := range entityLumpStrings {
        // skip empty entity strings
        if len(strings.TrimSpace(entityLumpString)) == 0 {
            continue
        }
        kvMap := make(map[string]string)
        linesScanner := bufio.NewScanner(strings.NewReader(entityLumpString))
        for linesScanner.Scan() {
            line := linesScanner.Text()
            if line == "{" || line == "}" || len(line) == 0 {
                continue
            }
            prop := strings.Split(line, "\" \"")
            key := strings.Trim(prop[0], "\"")
            value := strings.Trim(prop[1], "\"")
            kvMap[key] = value
        }
        entities = append(entities, kvMap)
    }

    outfile, err := os.Create("out.js")
    if err != nil {
        log.Fatal(err)
    }
    defer outfile.Close()

    for _, entity := range entities {
        classname := entity["classname"]
        if !strings.HasPrefix(classname, "trigger_") {
            continue
        }
        entityOriginParts := strings.Split(entity["origin"], " ")
        entityOriginX, _ := strconv.Atoi(entityOriginParts[0])
        entityOriginY, _ := strconv.Atoi(entityOriginParts[1])
        entityOriginZ, _ := strconv.Atoi(entityOriginParts[2])
        modelIndexStr := strings.Trim(entity["model"], "*")
        modelIndex, _ := strconv.Atoi(modelIndexStr)
        model := bspData.models[modelIndex]
        mesh := NewMeshFromBrushEntity()
        mesh.name = entity["targetname"]
        mesh.mins = model.Mins
        mesh.maxs = model.Maxs
        mesh.origin = model.Origin.Add(mgl32.Vec3{float32(entityOriginX), float32(entityOriginY), float32(entityOriginZ)})
        mesh.WriteToBabylon(outfile)
    }

    // all rendered faces go into a single mesh
    bspMesh := NewMesh()

    for _, face := range bspData.faces {
        stringIndex := int(bspData.texInfos[face.TexInfo].TexData)
        materialName, err := stringTable.FindString(stringIndex)
        if err != nil {
            log.Fatal(err)
        }
        if strings.HasPrefix(materialName, "TOOLS/") {
            continue
        }

        if face.DispInfo > -1 {
            //panic("don't want to handle displacements")
        } else {
           generateBspFace(&face, &bspData, bspMesh)
        }
    }

    // Write the faces to a Babylon mesh
    bspMesh.WriteToBabylon(outfile)

    // Get the planes used for collisions from the brush data, though
    for _, brush := range bspData.brushes {
        if (brush.Contents & bsp.MASK_PLAYERSOLID == 0) {
            continue
        }
        isPlayerClip := brush.Contents & bsp.CONTENTS_PLAYERCLIP != 0
        brushSides := bspData.brushSides[brush.FirstSide:(brush.FirstSide + brush.NumSides)]

        // check if this brush is 100% textured by tools
        // if so, then just skip it
        // note: player clip brushes get compiled into nodraw textured faces
        if (!isPlayerClip) {
            allToolTextures := true
            for _, side := range brushSides {
                stringIndex := int(bspData.texInfos[side.TexInfo].TexData)
                materialName, err := stringTable.FindString(stringIndex)
                if err != nil {
                    log.Fatal(err)
                }
                isAnyToolsTexture := strings.HasPrefix(materialName, "TOOLS/")
                if isPlayerClip {
                    fmt.Println(materialName)
                }
                if !isAnyToolsTexture || isPlayerClip {
                    allToolTextures = false
                    break
                }
            }
            if (allToolTextures) {
                continue
            }
        }

        planes := make([]*plane.Plane, len(brushSides))
        verticesBySide := make([][]mgl32.Vec3, len(brushSides))

        for sideIdx, side := range brushSides {
            verticesBySide[sideIdx] = make([]mgl32.Vec3, 0)
            planes[sideIdx] = &bspData.planes[side.PlaneNum]
        }

        meshFromPlanes := NewMeshFromPlanes()
        meshFromPlanes.planes = planes

        if (!brushCollisionsOnly) {
            // The idea here is to convert the planes into a mesh by:
            // - Getting the vertexes by computing the intersections of all combinations of 3 planes
            // - Converting the vertexes into 2d coorindates on a plane-local coordinate system
            // - Triangulating the 2d face into triangles via the triangolatte library
            // - Creating a mesh using the triangles to determine the indices into the vertex/normal array
            // This is probably dumb, though, since the data for the faces already exist in the bsp data.
            // The nice thing is that the brush data is split up into easier-to-understand shapes
            // so I thought this might be more managable even if its overly complicated/unnecessary.
            //
            // Ultimately, it needs more work to be functional, though, because it can't handle complex shapes
            // and float math leads to some weirdness.

            planeCombinations := combin.Combinations(len(planes), 3)
            for _, planeCombination := range planeCombinations {
                a, b, c := planeCombination[0], planeCombination[1], planeCombination[2]
                intersection, err := getIntersection(planes[a], planes[b], planes[c])
                if (err == nil) {
                    verticesBySide[a] = append(verticesBySide[a], intersection)
                    verticesBySide[b] = append(verticesBySide[b], intersection)
                    verticesBySide[c] = append(verticesBySide[c], intersection)
                }
            }

            vertexOffset := uint(0)

            for sideIdx, sideVertices := range verticesBySide {
                sidePlane := planes[sideIdx]
                xAxis, yAxis := getAxesFromPlane(sidePlane)
                points := make([]triangolatte.Point, len(sideVertices))
                pointsToVertexMap := make(map[triangolatte.Point]uint)

                for vertIdx, vertex := range sideVertices {
                    // add vert and normal to mesh in advance since we know
                    // they will all be needed
                    meshFromPlanes.AddVertexAndNormal(vertex, sidePlane.Normal)
                    // dummy uv
                    meshFromPlanes.AddUV(0, 1)

                    point := vertexToPlaneCoords(vertex, xAxis, yAxis)
                    points[vertIdx] = point
                    pointsToVertexMap[point] = uint(vertIdx)
                }

                // triangolatte requires vertices to be sorted in clock-wise order
                SortVertices(points)

                triangles, err := triangolatte.Polygon(points)
                if (len(triangles) == 0 || err != nil) {
                    fmt.Println(triangles)
                    fmt.Println(points)
                    panic(fmt.Sprintf("Failed to triangulate polygon: %s", err))
                }

                faceVertexIndexes := make([]uint, len(triangles)/2)
                for i := 0; i < len(triangles); i += 2 {
                    triangleX := triangles[i]
                    triangleY := triangles[i+1]
                    trianglePoint := triangolatte.Point{triangleX, triangleY}
                    vertexIdx := pointsToVertexMap[trianglePoint]
                    faceVertexIndexes[i/2] = vertexIdx
                    meshFromPlanes.AddIndices(vertexOffset + vertexIdx)
                }

                vertexOffset += uint(len(sideVertices))
            }
        }

        meshFromPlanes.WriteToBabylon(outfile, brushCollisionsOnly)
    }
}

func f32Abs(v float32) float32 {
    if v < 0 {
        return -v
    }
    return v
}

// Choose arbitrary x and y axes along the given plane
func getAxesFromPlane(plane3d *plane.Plane) (xAxis mgl32.Vec3, yAxis mgl32.Vec3) {
    // try projecting a regular x axis onto the plane
    xAxis = mgl32.Vec3{1,0,0}
    // handle normals that are parallel to our chosen x axis
    if (mgl32.FloatEqual(f32Abs(plane3d.Normal.X()), 1)) {
        xAxis = mgl32.Vec3{0,0,1}
    }
    xAxis = xAxis.Sub(plane3d.Normal.Mul(xAxis.Dot(plane3d.Normal)))
    xAxis = xAxis.Normalize()

    yAxis = plane3d.Normal.Cross(xAxis)
    return xAxis, yAxis
}

// Convert a point along a plane to the 2d coordinate system defined by the x and y axes
func vertexToPlaneCoords(vertex mgl32.Vec3, xAxis mgl32.Vec3, yAxis mgl32.Vec3) (triangolatte.Point) {
    x := vertex.Dot(xAxis)
    y := vertex.Dot(yAxis)
    return triangolatte.Point{float64(x), float64(y)}
}

func getIntersection(plane1 *plane.Plane, plane2 *plane.Plane, plane3 *plane.Plane) (intersection mgl32.Vec3, err error) {
    matrix := mgl32.Mat3FromCols(plane1.Normal, plane2.Normal, plane3.Normal)
    det := matrix.Det()

    if (det == 0) {
        return mgl32.Vec3{0,0,0}, errors.New("no intersection")
    }

    a := plane2.Normal.Cross(plane3.Normal).Mul(-plane1.Distance)
    b := plane3.Normal.Cross(plane1.Normal).Mul(-plane2.Distance)
    c := plane1.Normal.Cross(plane2.Normal).Mul(-plane3.Distance)
    return (a.Add(b).Add(c)).Mul(1/det), nil // (a+b+c)/det
}

func FindCentroid(points []triangolatte.Point) (triangolatte.Point) {
    x := 0.0
    y := 0.0
    for _, point := range points {
        x += point.X
        y += point.Y
    }
    center := triangolatte.Point{0,0}
    center.X = x / float64(len(points))
    center.Y = y / float64(len(points))
    return center
}

func SortVertices(points []triangolatte.Point) ([]triangolatte.Point) {
    center := FindCentroid(points)
    sort.Slice(points, func(i, j int) bool {
        a := points[i]
        b := points[j]
        a1 := math.Mod(ToDegrees(math.Atan2(a.X - center.X, a.Y - center.Y)) + 360, 360)
        a2 := math.Mod(ToDegrees(math.Atan2(b.X - center.X, b.Y - center.Y)) + 360, 360)
        return a1 > a2
    })
    return points
}

func ToDegrees(radians float64) float64 {
    return radians * (180.0 / math.Pi)
}

type MeshFromBrushEntity struct {
    mins            mgl32.Vec3
    maxs            mgl32.Vec3
    origin          mgl32.Vec3
    name            string
}

func NewMeshFromBrushEntity() *MeshFromBrushEntity {
    return &MeshFromBrushEntity{}
}

func (mesh *MeshFromBrushEntity) WriteToBabylon(writer io.StringWriter) {
    width := f32Abs(mesh.maxs.X() - mesh.mins.X())
    depth := f32Abs(mesh.maxs.Y() - mesh.mins.Y())
    height := f32Abs(mesh.maxs.Z() - mesh.mins.Z())

    meshDecl := fmt.Sprintf("var mesh = new BABYLON.MeshBuilder.CreateBox(\"%s\", {width: %f, depth: %f, height: %f}, scene);\n", mesh.name, width, depth, height);
    writer.WriteString(meshDecl)

    writer.WriteString("mesh.bspPlanes = [];\n")
    meshPos := fmt.Sprintf("mesh.position = new BABYLON.Vector3(%f, %f, %f);\n", mesh.origin.X(), mesh.origin.Z(), mesh.origin.Y())
    writer.WriteString(meshPos)
    // TODO: remove this, just for debugging purposes
    writer.WriteString("mesh.visibility = 0.25;\n")
    writer.WriteString("scene.triggers.push(mesh);\n")
}

type MeshFromPlanes struct {
    vertices            []float32
    normals             []float32
    uvs                 []float32
    indices             []uint
    planes              []*plane.Plane
}

func NewMeshFromPlanes() *MeshFromPlanes {
    return &MeshFromPlanes{}
}

func (mesh *MeshFromPlanes) AddVertexAndNormal(vertex mgl32.Vec3, normal mgl32.Vec3) {
    mesh.vertices = append(mesh.vertices, vertex.X(), vertex.Y(), vertex.Z())
    mesh.normals = append(mesh.normals, normal.X(), normal.Y(), normal.Z())
}

func (mesh *MeshFromPlanes) AddUV(uv ...float32) {
    mesh.uvs = append(mesh.uvs, uv...)
}

func (mesh *MeshFromPlanes) AddIndices(indices ...uint) {
    mesh.indices = append(mesh.indices, indices...)
}

func (mesh *MeshFromPlanes) ToBabylonCoordinates() *MeshFromPlanes {
    converted := NewMeshFromPlanes()
    // swap y and z, since y is up/down in babylon
    converted.vertices = make([]float32, len(mesh.vertices))
    converted.normals = make([]float32, len(mesh.normals))
    for i := uint(0); i < uint(len(mesh.vertices)); i += 3 {
        converted.vertices[i] = -mesh.vertices[i]
        converted.vertices[i+1] = -mesh.vertices[i+2]
        converted.vertices[i+2] = -mesh.vertices[i+1]
        converted.normals[i] = mesh.normals[i]
        converted.normals[i+1] = mesh.normals[i+2]
        converted.normals[i+2] = mesh.normals[i+1]
    }
    converted.uvs = mesh.uvs
    converted.indices = mesh.indices
    converted.planes = mesh.planes
    return converted
}

func (mesh *MeshFromPlanes) WriteToBabylon(writer io.StringWriter, collisionsOnly bool) {
    converted := mesh.ToBabylonCoordinates()

    bspId += 1

    meshDecl := fmt.Sprintf("var mesh = new BABYLON.Mesh(\"bsp%d\", scene);\n", bspId);
    writer.WriteString(meshDecl)

    if (!collisionsOnly) {
        writer.WriteString("var vertexData = new BABYLON.VertexData();\n")

        writer.WriteString("vertexData.indices = ")
        indicesJson, err := json.Marshal(converted.indices)
        if (err != nil) {
            log.Fatal(err)
        }
        writer.WriteString(string(indicesJson))
        writer.WriteString(";\n")

        writer.WriteString("vertexData.positions = ")
        verticesJson, err := json.Marshal(converted.vertices)
        if (err != nil) {
            log.Fatal(err)
        }
        writer.WriteString(string(verticesJson))
        writer.WriteString(";\n")

        writer.WriteString("vertexData.normals = ")
        normalsJson, err := json.Marshal(converted.normals)
        if (err != nil) {
            log.Fatal(err)
        }
        writer.WriteString(string(normalsJson))
        writer.WriteString(";\n")

        writer.WriteString("vertexData.uvs = ")
        uvsJson, err := json.Marshal(converted.uvs)
        if (err != nil) {
            log.Fatal(err)
        }
        writer.WriteString(string(uvsJson))
        writer.WriteString(";\n")

        writer.WriteString("vertexData.applyToMesh(mesh, false);\n")
    }

    writer.WriteString("mesh.bspPlanes = [];\n")
    for _, sidePlane := range mesh.planes {
        dist := strconv.FormatFloat(float64(sidePlane.Distance), 'f', -1, 32)
        x := strconv.FormatFloat(float64(sidePlane.Normal[0]), 'f', -1, 32)
        // y is up/down in babylon
        y := strconv.FormatFloat(float64(sidePlane.Normal[2]), 'f', -1, 32)
        z := strconv.FormatFloat(float64(sidePlane.Normal[1]), 'f', -1, 32)
        s := fmt.Sprintf("mesh.bspPlanes.push({dist: %s, normal: new BABYLON.Vector3(%s, %s, %s)});\n", dist, x, y, z)
        writer.WriteString(s)
    }
}

// generateBspFace Create primitives from face data in the bsp
func generateBspFace(f *face.Face, bspData *BSPData, bspMesh IMesh) {
    //offset := int32(len(bspMesh.Vertices())) / 3
    length := int32(0)

    planeNormal := bspData.planes[f.Planenum].Normal
    // All surfedges associated with this face
    // surfEdges are basically indices into the edges lump
    faceSurfEdges := bspData.surfEdges[f.FirstEdge:(f.FirstEdge + int32(f.NumEdges))]
    rootIndex := uint16(0)
    for idx, surfEdge := range faceSurfEdges {
        edge := bspData.edges[int(math.Abs(float64(surfEdge)))]
        e1 := 0
        e2 := 1
        if surfEdge < 0 {
            e1 = 1
            e2 = 0
        }
        //Capture root indice
        if idx == 0 {
            rootIndex = edge[e1]
        } else {
            // Just create a triangle for every edge now
            bspMesh.AddVertex(bspData.vertexes[rootIndex].X(), bspData.vertexes[rootIndex].Y(), bspData.vertexes[rootIndex].Z())
            bspMesh.AddNormal(planeNormal.X(), planeNormal.Y(), planeNormal.Z())

            bspMesh.AddVertex(bspData.vertexes[edge[e1]].X(), bspData.vertexes[edge[e1]].Y(), bspData.vertexes[edge[e1]].Z())
            bspMesh.AddNormal(planeNormal.X(), planeNormal.Y(), planeNormal.Z())

            bspMesh.AddVertex(bspData.vertexes[edge[e2]].X(), bspData.vertexes[edge[e2]].Y(), bspData.vertexes[edge[e2]].Z())
            bspMesh.AddNormal(planeNormal.X(), planeNormal.Y(), planeNormal.Z())

            length += 3 // num verts (3 b/c face triangles
        }
    }
}

// IMesh Generic Mesh interface
// Most renderable objects should implement this, but there
// are probably many custom cases that may not
type IMesh interface {
    // AddVertex
    AddVertex(...float32)
    // AddNormal
    AddNormal(...float32)
    // AddUV
    AddUV(...float32)

    // Vertices
    Vertices() []float32
    // Normals
    Normals() []float32
    // UVs
    UVs() []float32
}

// Mesh
type Mesh struct {
    vertices            []float32
    normals             []float32
    uvs                 []float32
}

// AddVertex
func (mesh *Mesh) AddVertex(vertex ...float32) {
    mesh.vertices = append(mesh.vertices, vertex...)
}

// AddNormal
func (mesh *Mesh) AddNormal(normal ...float32) {
    mesh.normals = append(mesh.normals, normal...)
}

// AddUV
func (mesh *Mesh) AddUV(uv ...float32) {
    mesh.uvs = append(mesh.uvs, uv...)
}

// Vertices
func (mesh *Mesh) Vertices() []float32 {
    return mesh.vertices
}

// Normals
func (mesh *Mesh) Normals() []float32 {
    return mesh.normals
}

// UVs
func (mesh *Mesh) UVs() []float32 {
    return mesh.uvs
}

func (mesh *Mesh) ToBabylonCoordinates() *Mesh {
    converted := NewMesh()
    // swap y and z, since y is up/down in babylon
    converted.vertices = make([]float32, len(mesh.vertices))
    converted.normals = make([]float32, len(mesh.normals))
    for i := uint(0); i < uint(len(mesh.vertices)); i += 3 {
        converted.vertices[i] = mesh.vertices[i]
        converted.vertices[i+1] = mesh.vertices[i+2]
        converted.vertices[i+2] = mesh.vertices[i+1]
        converted.normals[i] = mesh.normals[i]
        converted.normals[i+1] = mesh.normals[i+2]
        converted.normals[i+2] = mesh.normals[i+1]
    }
    converted.uvs = mesh.uvs
    return converted
}

func (mesh *Mesh) WriteToBabylon(writer io.StringWriter) {
    converted := mesh.ToBabylonCoordinates()

    bspId += 1

    numVertices := len(mesh.vertices)/3
    numIndices := numVertices
    indices := make([]uint, numIndices)
    uvs := make([]float32, numVertices * 2)
    for i := uint(0); i < uint(numIndices); i += 1 {
        indices[i] = uint(numIndices-1)-i
        uvs[i*2] = 0
        uvs[i*2+1] = 1
    }

    meshDecl := fmt.Sprintf("var mesh = new BABYLON.Mesh(\"bsp%d\", scene);\n", bspId);
    writer.WriteString(meshDecl)
    writer.WriteString("var vertexData = new BABYLON.VertexData();\n")

    writer.WriteString("vertexData.indices = ")
    indicesJson, err := json.Marshal(indices)
    if (err != nil) {
        log.Fatal(err)
    }
    writer.WriteString(string(indicesJson))
    writer.WriteString(";\n")

    writer.WriteString("vertexData.positions = ")
    verticesJson, err := json.Marshal(converted.vertices)
    if (err != nil) {
        log.Fatal(err)
    }
    writer.WriteString(string(verticesJson))
    writer.WriteString(";\n")

    writer.WriteString("vertexData.normals = ")
    normalsJson, err := json.Marshal(converted.normals)
    if (err != nil) {
        log.Fatal(err)
    }
    writer.WriteString(string(normalsJson))
    writer.WriteString(";\n")

    writer.WriteString("vertexData.uvs = ")
    uvsJson, err := json.Marshal(uvs)
    if (err != nil) {
        log.Fatal(err)
    }
    writer.WriteString(string(uvsJson))
    writer.WriteString(";\n")

    writer.WriteString("vertexData.applyToMesh(mesh, false);\n")
}

// NewMesh
func NewMesh() *Mesh {
    return &Mesh{}
}
