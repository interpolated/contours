import { useEffect, useState } from "react";
import { useGiraffeState } from "@gi-nx/iframe-sdk-react";
import { rpc } from "@gi-nx/iframe-sdk";
import proj4 from "proj4";
import * as THREE from "three";
import * as turf from "@turf/turf";
import { GLTFExporter } from 'three-stdlib';

// Helper to create projection
const getReferenceGeoProject = (ref: [number, number]) =>
  proj4(`+proj=lcc +lat_0=${ref[1]} +lat_1=${ref[1] + 0.1} +lat_2=${ref[1] - 0.1} +lon_0=${ref[0]} +units=m +no_defs +ellps=WGS84`);

function App() {
  const origin = useGiraffeState("projectOrigin");
  const baked = useGiraffeState("bakedSections");
  const [status, setStatus] = useState("");

  const runMesh = async () => {
    if (!origin || !baked?.features) return;

    setStatus("Sampling points...");
    const { forward } = getReferenceGeoProject(origin);

    const contours = baked.features.filter(f => f.properties?.contour);
    const sampledPoints = [];

    for (const feature of contours) {
      const line = turf.lineString(feature.geometry.coordinates);
      const length = turf.length(line, { units: "meters" });
      for (let i = 0; i <= length; i += 0.005) {
        const pt = turf.along(line, i, { units: "kilometers" });
        const [lng, lat] = pt.geometry.coordinates;
        const projected = forward([lng, lat]);
        sampledPoints.push({
          geometry: { x: projected[0], y: projected[1] },
          properties: { elevation: feature.properties.elevation ?? 0 }
        });
      }
    }

    if (sampledPoints.length === 0) {
      setStatus("No points sampled.");
      return;
    }

    setStatus("Creating mesh...");
    const meshGroup = createTriangularMesh(sampledPoints);

    const scene = new THREE.Scene();
    meshGroup.rotateX(-Math.PI / 2);
    scene.add(meshGroup);

    const exporter = new GLTFExporter();
    setStatus("Exporting GLB...");
    exporter.parse(
      scene,
      async function (result) {
        if (result instanceof ArrayBuffer) {
          const dataUrl = await arrayBufferToDataUrl(result);
          await saveGltfAndAddLayer(dataUrl);
          setStatus("Export complete.");
        } else {
          console.error("Expected ArrayBuffer for GLB export");
          setStatus("Export failed.");
        }
      },
      (error) => {
        console.error("Export error", error);
        setStatus("Export failed.");
      },
      { binary: true }
    );
  };

  async function saveGltfAndAddLayer(modelUri: string) {
    const style = {
      id: "landscape",
      type: "model",
      slot: "middle",
      source: {
        type: "geojson",
        data: {
          type: "Feature",
          properties: { "model-uri": modelUri },
          geometry: {
            type: "Point",
            coordinates: [origin[0], origin[1]]
          }
        }
      },
      minzoom: 15,
      layout: { "model-id": ["get", "model-uri"] },
      paint: {
        "model-opacity": 1,
        "model-rotation": [0, 0, 0],
        "model-scale": [1, 1, 1],
        "model-color-mix-intensity": 0,
        "model-cast-shadows": true,
        "model-receive-shadows": true,
        "model-emissive-strength": 0,
        "model-height-based-emissive-strength-multiplier": [0, 10, 1, 100, 1],
        "model-ambient-occlusion-intensity": 1
      }
    };

    await rpc.invoke("addTempLayer", ["mesh", style]);
  }

  async function arrayBufferToDataUrl(buffer: ArrayBuffer): Promise<string> {
    return new Promise((resolve, reject) => {
      const blob = new Blob([buffer], { type: "model/gltf-binary" });
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result as string);
      reader.onerror = reject;
      reader.readAsDataURL(blob);
    });
  }

  return (
    <div style={{ padding: "20px" }}>
      <button onClick={runMesh}>Export Triangular Mesh</button>
      <p>{status}</p>
    </div>
  );
}

export default App;