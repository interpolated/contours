import { useCallback, useEffect, useRef, useState } from "react";
import { useGiraffeState } from "@gi-nx/iframe-sdk-react";
import { rpc } from "@gi-nx/iframe-sdk";
import proj4 from "proj4";
import * as THREE from "three";
import * as turf from "@turf/turf";
import { GLTFExporter } from 'three-stdlib';


import mapboxgl from "mapbox-gl";

// Helper to create projection
const getReferenceGeoProject = (ref) =>
  proj4(`+proj=lcc +lat_0=${ref[1]} +lat_1=${ref[1] + 0.1} +lat_2=${ref[1] - 0.1} +lon_0=${ref[0]} +units=m +no_defs +ellps=WGS84`);


mapboxgl.accessToken = 'pk.eyJ1IjoiZ2FydGhkYmVldGxlIiwiYSI6ImNpcHl5emhrdjB5YmxoY25yczF6MHhhc2IifQ.2Ld30uLqcffVv-RUAWk_qQ'; // Replace with your token


function App() {
    const origin = useGiraffeState("projectOrigin");
    const baked = useGiraffeState("bakedSections");
    const mapContainer = useRef<HTMLDivElement>(null);
    const mapRef = useRef<mapboxgl.Map | null>(null);
    const [modelAdded, setModelAdded] = useState(false);
  
    useEffect(() => {
      if (mapContainer.current && !mapRef.current && origin) {
        const map = new mapboxgl.Map({
          container: mapContainer.current,
          style: "mapbox://styles/mapbox/streets-v12",
          center: origin,
          zoom: 16,
          pitch: 60,
          bearing: -17.6,
          antialias: true
        });
  
        map.on("load", () => {
          mapRef.current = map;
        });
      }
    }, [origin]);
  
    const runMesh = async () => {
      const { forward } = getReferenceGeoProject(origin);
  
      const contours = baked?.features?.filter(f => f.properties?.contour) || [];
      const sampledPoints = [];
      const originalPoints = [];
  
      for (const feature of contours) {
        const line = turf.lineString(feature.geometry.coordinates);
        const length = turf.length(line, { units: "meters" });
        for (let i = 0; i <= length; i += 0.005) {
          const pt = turf.along(line, i, { units: "kilometers" });
          const [lng, lat] = pt.geometry.coordinates;
          originalPoints.push([lng, lat]);
          const projected = forward([lng, lat]);
          sampledPoints.push({
            geometry: { x: projected[0], y: projected[1] },
            properties: { elevation: feature.properties.elevation ?? 0 }
          });
        }
      }
  
      if (sampledPoints.length === 0) {
        console.warn("No points sampled");
        return;
      }
  
      if (mapRef.current && originalPoints.length > 0) {
        const bbox = originalPoints.reduce(
          (bounds, coord) => bounds.extend(coord as [number, number]),
          new mapboxgl.LngLatBounds(originalPoints[0] as [number, number], originalPoints[0] as [number, number])
        );
        mapRef.current.fitBounds(bbox, { padding: 50, duration: 1000 });
      }
  
      const rawMeshGroup = createTriangularMesh(sampledPoints);
  
      let meshGroup;
      if (rawMeshGroup instanceof THREE.Object3D) {
        meshGroup = rawMeshGroup;
      } else {
        const geometry = rawMeshGroup;
        const material = new THREE.MeshStandardMaterial({ color: 0xcccccc, side: THREE.DoubleSide });
        meshGroup = new THREE.Mesh(geometry, material);
      }
      const scene = new THREE.Scene();

      meshGroup.rotateX(-Math.PI / 2);

      scene.add(meshGroup);
  
      const exporter = new GLTFExporter();
      exporter.parse(
        scene,
        async function (result) {
          if (result instanceof ArrayBuffer) {
            const dataUrl = await arrayBufferToDataUrl(result);
            await saveGltfAndAddLayer(dataUrl);
          } else {
            console.error("Expected ArrayBuffer for GLB export");
          }
        },
        function (error) {
          console.error("An error happened during export", error);
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
            properties: {
              "model-uri": modelUri
            },
            geometry: {
              type: "Point",
              coordinates: [origin[0], origin[1]]
            }
          }
        },
        minzoom: 15,
        layout: {
          "model-id": ["get", "model-uri"]
        },
        paint: {
          "model-opacity": 1,
          "model-rotation": [0 , 0, 0],
          "model-scale": [1, 1, 1],
          "model-color-mix-intensity": 0,
          "model-cast-shadows": true,
            "model-receive-shadows": true,
          "model-emissive-strength": 0,
          "model-height-based-emissive-strength-multiplier":[0,10,1,100,1],
          "model-ambient-occlusion-intensity": 1
        }
      };
  
      console.log("Adding model to Mapbox map");
  
      if (mapRef.current && !modelAdded) {
        const map = mapRef.current;
      
        if (map.getLayer("landscape")) {
          map.removeLayer("landscape");
        }
    
        if (map.getSource("landscape-source3")) {
          map.removeSource("landscape-source3");
        }
      
        map.addSource("landscape-source3", style.source);
      
        map.addLayer({
          ...style,
          source: "landscape-source3"
        });
      
      }

      await rpc.invoke('addTempLayer', ['mesh', style]);
    }
  
    async function arrayBufferToDataUrl(buffer: ArrayBuffer): Promise<string> {
      return new Promise((resolve, reject) => {
        const blob = new Blob([buffer], { type: "model/gltf-binary" });
        const reader = new FileReader();
        reader.onload = () => {
          resolve(reader.result as string);
        };
        reader.onerror = reject;
        reader.readAsDataURL(blob);
      });
    }
  
    return (
      <div style={{ width: "100%", height: "100vh", position: "relative" }}>
        <div ref={mapContainer} style={{ width: "100%", height: "100%" }} />
        <button
          onClick={runMesh}
          style={{
            position: "absolute",
            top: 20,
            left: 20,
            zIndex: 1,
            padding: "10px",
            backgroundColor: "white",
            borderRadius: "5px"
          }}
        >
          Export and Add Model
        </button>
      </div>
    );
  }
  

// Helper functions and constants
const e = 11102230246251565e-32,
      n = 134217729,
      o = (3 + 8 * e) * e;

function r(t, e, n, o, r) {
    let f, i, u, s, c = e[0], a = o[0], l = 0, d = 0;
    a > c == a > -c ? (f = c, c = e[++l]) : (f = a, a = o[++d]);
    let p = 0;
    if (l < t && d < n) {
        for (a > c == a > -c ? (i = c + f, u = f - (i - c), c = e[++l]) : (i = a + f, u = f - (i - a), a = o[++d]), f = i, 0 !== u && (r[p++] = u); l < t && d < n;) {
            a > c == a > -c ? (i = f + c, s = i - f, u = f - (i - s) + (c - s), c = e[++l]) : (i = f + a, s = i - f, u = f - (i - s) + (a - s), a = o[++d]), f = i, 0 !== u && (r[p++] = u);
        }
        for (; l < t;) i = f + c, s = i - f, u = f - (i - s) + (c - s), c = e[++l], f = i, 0 !== u && (r[p++] = u);
        for (; d < n;) i = f + a, s = i - f, u = f - (i - s) + (a - s), a = o[++d], f = i, 0 !== u && (r[p++] = u);
    }
    return 0 === f && 0 !== p || (r[p++] = f), p;
}

function f(t) {
    return new Float64Array(t);
}




// Define orient2d function
function orient2d(t, e, f, p, b, h) {
    const y = (e - h) * (f - b),
          x = (t - b) * (p - h),
          M = y - x,
          g = Math.abs(y + x);
    return Math.abs(M) >= 33306690738754716e-32 * g ? M : -function(t, e, f, p, b, h, y) {
        let x, M, g, m, T, j, w, A, F, k, q, v, z, B, C, D, E, G;
        const H = t - b, I = f - b, J = e - h, K = p - h;
        B = H * K, j = n * H, w = j - (j - H), A = H - w, j = n * K, F = j - (j - K), k = K - F;
        C = A * k - (B - w * F - A * F - w * k), D = J * I, j = n * J, w = j - (j - J), A = J - w, j = n * I, F = j - (j - I), k = I - F;
        E = A * k - (D - w * F - A * F - w * k), q = C - E, T = C - q, s[0] = C - (q + T) + (T - E), v = B + q, T = v - B, z = B - (v - T) + (q - T), q = z - D, T = z - q, s[1] = z - (q + T) + (T - D), G = v + q, T = G - v, s[2] = v - (G - T) + (q - T), s[3] = G;

        let L = function(t, e) {
            let n = e[0];
            for (let o = 1; o < t; o++) n += e[o];
            return n;
        }(4, s), N = i * y;
        if (L >= N || -L >= N) return L;
        if (T = t - H, x = t - (H + T) + (T - b), T = f - I, g = f - (I + T) + (T - b), T = e - J, M = e - (J + T) + (T - h), T = p - K, m = p - (K + T) + (T - h), 0 === x && 0 === M && 0 === g && 0 === m) return L;
        if (N = u * y + o * Math.abs(L), L += H * m + K * x - (J * g + I * M), L >= N || -L >= N) return L;
        B = x * K, j = n * x, w = j - (j - x), A = x - w, j = n * K, F = j - (j - K), k = K - F;
        C = A * k - (B - w * F - A * F - w * k), D = M * I, j = n * M, w = j - (j - M), A = M - w, j = n * I, F = j - (j - I), k = I - F;
        E = A * k - (D - w * F - A * F - w * k), q = C - E, T = C - q, d[0] = C - (q + T) + (T - E), v = B + q, T = v - B, z = B - (v - T) + (q - T), q = z - D, T = z - q, d[1] = z - (q + T) + (T - D), G = v + q, T = G - v, d[2] = v - (G - T) + (q - T), d[3] = G;

        const O = r(4, s, 4, d, c);
        B = H * m, j = n * H, w = j - (j - H), A = H - w, j = n * m, F = j - (j - m), k = m - F;
        C = A * k - (B - w * F - A * F - w * k), D = J * g, j = n * J, w = j - (j - J), A = J - w, j = n * g, F = j - (j - g), k = g - F;
        E = A * k - (D - w * F - A * F - w * k), q = C - E, T = C - q, d[0] = C - (q + T) + (T - E), v = B + q, T = v - B, z = B - (v - T) + (q - T), q = z - D, T = z - q, d[1] = z - (q + T) + (T - D), G = v + q, T = G - v, d[2] = v - (G - T) + (q - T), d[3] = G;

        const P = r(O, c, 4, d, a);
        return a[P - 1];
    }(t, e, f, p, b, h, y);
}


const EPSILON = Math.pow(2, -52);
const EDGE_STACK = new Uint32Array(512);



class Delaunator {

    static from(points, getX = defaultGetX, getY = defaultGetY) {
        const n = points.length;
        const coords = new Float64Array(n * 2);

        for (let i = 0; i < n; i++) {
            const p = points[i];
            coords[2 * i] = getX(p);
            coords[2 * i + 1] = getY(p);
        }

        return new Delaunator(coords);
    }

    constructor(coords) {
        const n = coords.length >> 1;
        if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.');

        this.coords = coords;

        // arrays that will store the triangulation graph
        const maxTriangles = Math.max(2 * n - 5, 0);
        this._triangles = new Uint32Array(maxTriangles * 3);
        this._halfedges = new Int32Array(maxTriangles * 3);

        // temporary arrays for tracking the edges of the advancing convex hull
        this._hashSize = Math.ceil(Math.sqrt(n));
        this._hullPrev = new Uint32Array(n); // edge to prev edge
        this._hullNext = new Uint32Array(n); // edge to next edge
        this._hullTri = new Uint32Array(n); // edge to adjacent triangle
        this._hullHash = new Int32Array(this._hashSize); // angular edge hash

        // temporary arrays for sorting points
        this._ids = new Uint32Array(n);
        this._dists = new Float64Array(n);

        this.update();
    }

    update() {
        const {coords, _hullPrev: hullPrev, _hullNext: hullNext, _hullTri: hullTri, _hullHash: hullHash} =  this;
        const n = coords.length >> 1;

        // populate an array of point indices; calculate input data bbox
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;

        for (let i = 0; i < n; i++) {
            const x = coords[2 * i];
            const y = coords[2 * i + 1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
            this._ids[i] = i;
        }
        const cx = (minX + maxX) / 2;
        const cy = (minY + maxY) / 2;

        let i0, i1, i2;

        // pick a seed point close to the center
        for (let i = 0, minDist = Infinity; i < n; i++) {
            const d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
            if (d < minDist) {
                i0 = i;
                minDist = d;
            }
        }
        const i0x = coords[2 * i0];
        const i0y = coords[2 * i0 + 1];

        // find the point closest to the seed
        for (let i = 0, minDist = Infinity; i < n; i++) {
            if (i === i0) continue;
            const d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
            if (d < minDist && d > 0) {
                i1 = i;
                minDist = d;
            }
        }
        let i1x = coords[2 * i1];
        let i1y = coords[2 * i1 + 1];

        let minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (let i = 0; i < n; i++) {
            if (i === i0 || i === i1) continue;
            const r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
            if (r < minRadius) {
                i2 = i;
                minRadius = r;
            }
        }
        let i2x = coords[2 * i2];
        let i2y = coords[2 * i2 + 1];

        if (minRadius === Infinity) {
            // order collinear points by dx (or dy if all x are identical)
            // and return the list as a hull
            for (let i = 0; i < n; i++) {
                this._dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
            }
            quicksort(this._ids, this._dists, 0, n - 1);
            const hull = new Uint32Array(n);
            let j = 0;
            for (let i = 0, d0 = -Infinity; i < n; i++) {
                const id = this._ids[i];
                const d = this._dists[id];
                if (d > d0) {
                    hull[j++] = id;
                    d0 = d;
                }
            }
            this.hull = hull.subarray(0, j);
            this.triangles = new Uint32Array(0);
            this.halfedges = new Uint32Array(0);
            return;
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0) {
            const i = i1;
            const x = i1x;
            const y = i1y;
            i1 = i2;
            i1x = i2x;
            i1y = i2y;
            i2 = i;
            i2x = x;
            i2y = y;
        }

        const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
        this._cx = center.x;
        this._cy = center.y;

        for (let i = 0; i < n; i++) {
            this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
        }

        // sort the points by distance from the seed triangle circumcenter
        quicksort(this._ids, this._dists, 0, n - 1);

        // set up the seed triangle as the starting hull
        this._hullStart = i0;
        let hullSize = 3;

        hullNext[i0] = hullPrev[i2] = i1;
        hullNext[i1] = hullPrev[i0] = i2;
        hullNext[i2] = hullPrev[i1] = i0;

        hullTri[i0] = 0;
        hullTri[i1] = 1;
        hullTri[i2] = 2;

        hullHash.fill(-1);
        hullHash[this._hashKey(i0x, i0y)] = i0;
        hullHash[this._hashKey(i1x, i1y)] = i1;
        hullHash[this._hashKey(i2x, i2y)] = i2;

        this.trianglesLen = 0;
        this._addTriangle(i0, i1, i2, -1, -1, -1);

        for (let k = 0, xp, yp; k < this._ids.length; k++) {
            const i = this._ids[k];
            const x = coords[2 * i];
            const y = coords[2 * i + 1];

            // skip near-duplicate points
            if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
            xp = x;
            yp = y;

            // skip seed triangle points
            if (i === i0 || i === i1 || i === i2) continue;

            // find a visible edge on the convex hull using edge hash
            let start = 0;
            for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
                start = hullHash[(key + j) % this._hashSize];
                if (start !== -1 && start !== hullNext[start]) break;
            }

            start = hullPrev[start];
            let e = start, q;
            while (q = hullNext[e], orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]) >= 0) {
                e = q;
                if (e === start) {
                    e = -1;
                    break;
                }
            }
            if (e === -1) continue; // likely a near-duplicate point; skip it

            // add the first triangle from the point
            let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i] = this._legalize(t + 2);
            hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize++;

            // walk forward through the hull, adding more triangles and flipping recursively
            let n = hullNext[e];
            while (q = hullNext[n], orient2d(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1]) < 0) {
                t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                hullTri[i] = this._legalize(t + 2);
                hullNext[n] = n; // mark as removed
                hullSize--;
                n = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e === start) {
                while (q = hullPrev[e], orient2d(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]) < 0) {
                    t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                    this._legalize(t + 2);
                    hullTri[q] = t;
                    hullNext[e] = e; // mark as removed
                    hullSize--;
                    e = q;
                }
            }

            // update the hull indices
            this._hullStart = hullPrev[i] = e;
            hullNext[e] = hullPrev[n] = i;
            hullNext[i] = n;

            // save the two new edges in the hash table
            hullHash[this._hashKey(x, y)] = i;
            hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
        }

        this.hull = new Uint32Array(hullSize);
        for (let i = 0, e = this._hullStart; i < hullSize; i++) {
            this.hull[i] = e;
            e = hullNext[e];
        }

        // trim typed triangle mesh arrays
        this.triangles = this._triangles.subarray(0, this.trianglesLen);
        this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
    }

    _hashKey(x, y) {
        return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
    }

    _legalize(a) {
        const {_triangles: triangles, _halfedges: halfedges, coords} = this;

        let i = 0;
        let ar = 0;

        // recursion eliminated with a fixed-size stack
        while (true) {
            const b = halfedges[a];

            /* if the pair of triangles doesn't satisfy the Delaunay condition
             * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
             * then do the same check/flip recursively for the new pair of triangles
             *
             *           pl                    pl
             *          /||\                  /  \
             *       al/ || \bl            al/    \a
             *        /  ||  \              /      \
             *       /  a||b  \    flip    /___ar___\
             *     p0\   ||   /p1   =>   p0\---bl---/p1
             *        \  ||  /              \      /
             *       ar\ || /br             b\    /br
             *          \||/                  \  /
             *           pr                    pr
             */
            const a0 = a - a % 3;
            ar = a0 + (a + 2) % 3;

            if (b === -1) { // convex hull edge
                if (i === 0) break;
                a = EDGE_STACK[--i];
                continue;
            }

            const b0 = b - b % 3;
            const al = a0 + (a + 1) % 3;
            const bl = b0 + (b + 2) % 3;

            const p0 = triangles[ar];
            const pr = triangles[a];
            const pl = triangles[al];
            const p1 = triangles[bl];

            const illegal = inCircle(
                coords[2 * p0], coords[2 * p0 + 1],
                coords[2 * pr], coords[2 * pr + 1],
                coords[2 * pl], coords[2 * pl + 1],
                coords[2 * p1], coords[2 * p1 + 1]);

            if (illegal) {
                triangles[a] = p1;
                triangles[b] = p0;

                const hbl = halfedges[bl];

                // edge swapped on the other side of the hull (rare); fix the halfedge reference
                if (hbl === -1) {
                    let e = this._hullStart;
                    do {
                        if (this._hullTri[e] === bl) {
                            this._hullTri[e] = a;
                            break;
                        }
                        e = this._hullPrev[e];
                    } while (e !== this._hullStart);
                }
                this._link(a, hbl);
                this._link(b, halfedges[ar]);
                this._link(ar, bl);

                const br = b0 + (b + 1) % 3;

                // don't worry about hitting the cap: it can only happen on extremely degenerate input
                if (i < EDGE_STACK.length) {
                    EDGE_STACK[i++] = br;
                }
            } else {
                if (i === 0) break;
                a = EDGE_STACK[--i];
            }
        }

        return ar;
    }

    _link(a, b) {
        this._halfedges[a] = b;
        if (b !== -1) this._halfedges[b] = a;
    }

    // add a new triangle given vertex indices and adjacent half-edge ids
    _addTriangle(i0, i1, i2, a, b, c) {
        const t = this.trianglesLen;

        this._triangles[t] = i0;
        this._triangles[t + 1] = i1;
        this._triangles[t + 2] = i2;

        this._link(t, a);
        this._link(t + 1, b);
        this._link(t + 2, c);

        this.trianglesLen += 3;

        return t;
    }
}

// monotonically increases with real angle, but doesn't need expensive trigonometry
function pseudoAngle(dx, dy) {
    const p = dx / (Math.abs(dx) + Math.abs(dy));
    return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
}

function dist(ax, ay, bx, by) {
    const dx = ax - bx;
    const dy = ay - by;
    return dx * dx + dy * dy;
}

function inCircle(ax, ay, bx, by, cx, cy, px, py) {
    const dx = ax - px;
    const dy = ay - py;
    const ex = bx - px;
    const ey = by - py;
    const fx = cx - px;
    const fy = cy - py;

    const ap = dx * dx + dy * dy;
    const bp = ex * ex + ey * ey;
    const cp = fx * fx + fy * fy;

    return dx * (ey * cp - bp * fy) -
           dy * (ex * cp - bp * fx) +
           ap * (ex * fy - ey * fx) < 0;
}

function circumradius(ax, ay, bx, by, cx, cy) {
    const dx = bx - ax;
    const dy = by - ay;
    const ex = cx - ax;
    const ey = cy - ay;

    const bl = dx * dx + dy * dy;
    const cl = ex * ex + ey * ey;
    const d = 0.5 / (dx * ey - dy * ex);

    const x = (ey * bl - dy * cl) * d;
    const y = (dx * cl - ex * bl) * d;

    return x * x + y * y;
}

function circumcenter(ax, ay, bx, by, cx, cy) {
    const dx = bx - ax;
    const dy = by - ay;
    const ex = cx - ax;
    const ey = cy - ay;

    const bl = dx * dx + dy * dy;
    const cl = ex * ex + ey * ey;
    const d = 0.5 / (dx * ey - dy * ex);

    const x = ax + (ey * bl - dy * cl) * d;
    const y = ay + (dx * cl - ex * bl) * d;

    return {x, y};
}

function quicksort(ids, dists, left, right) {
    if (right - left <= 20) {
        for (let i = left + 1; i <= right; i++) {
            const temp = ids[i];
            const tempDist = dists[temp];
            let j = i - 1;
            while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
            ids[j + 1] = temp;
        }
    } else {
        const median = (left + right) >> 1;
        let i = left + 1;
        let j = right;
        swap(ids, median, i);
        if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
        if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
        if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

        const temp = ids[i];
        const tempDist = dists[temp];
        while (true) {
            do i++; while (dists[ids[i]] < tempDist);
            do j--; while (dists[ids[j]] > tempDist);
            if (j < i) break;
            swap(ids, i, j);
        }
        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left) {
            quicksort(ids, dists, i, right);
            quicksort(ids, dists, left, j - 1);
        } else {
            quicksort(ids, dists, left, j - 1);
            quicksort(ids, dists, i, right);
        }
    }
}

function swap(arr, i, j) {
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function defaultGetX(p) {
    return p[0];
}
function defaultGetY(p) {
    return p[1];
}


// Function to convert your data to Delaunay-compatible format
function convertToDelaunayFormat(dataArray) {
  return dataArray.map(item => [item.geometry.x, item.geometry.y]);
}

// Function to convert your data to Delaunay-compatible format
function convertToDelaunayFormat3D(dataArray) {
  return dataArray.map(item => [item.geometry.x, item.geometry.y,item.properties.elevation]);
}

/////////////////////////working function

function createTriangularMesh(dataArray) {
    const points = convertToDelaunayFormat(dataArray);
    const points3d = convertToDelaunayFormat3D(dataArray);
  
    const delaunay = new Delaunator(points.flat());
    const triangles = delaunay.triangles;
  
    const positions = [];
    const indices = [];
    const colors = [];
    let baseIndex = 0;
  
    // First, find min and max z for normalization
    let minZ = Infinity;
    let maxZ = -Infinity;
    for (const p of points3d) {
      if (p[2] < minZ) minZ = p[2];
      if (p[2] > maxZ) maxZ = p[2];
    }
  
    const darkOrange = new THREE.Color(0xabbaad); // Dark orange
    const lightYellow = new THREE.Color(0xffffe0); // Light yellow
  
    for (let i = 0; i < triangles.length; i += 3) {
      const [a, b, c] = [triangles[i], triangles[i + 1], triangles[i + 2]];
  
      const p1 = new THREE.Vector3(...points3d[a]);
      const p2 = new THREE.Vector3(...points3d[b]);
      const p3 = new THREE.Vector3(...points3d[c]);
  
      const extrusionHeight = 0.3;
      const p1Extruded = new THREE.Vector3(p1.x, p1.y + extrusionHeight, p1.z);
      const p2Extruded = new THREE.Vector3(p2.x, p2.y + extrusionHeight, p2.z);
      const p3Extruded = new THREE.Vector3(p3.x, p3.y + extrusionHeight, p3.z);
  
      const vertices = [p1, p2, p3, p1Extruded, p2Extruded, p3Extruded];
  
      for (const v of vertices) {
        positions.push(v.x, v.y, v.z);
  
        const normalizedZ = (v.z - minZ) / (maxZ - minZ);
        const color = darkOrange.clone().lerp(lightYellow, normalizedZ);
        colors.push(color.r, color.g, color.b);
      }
  
      indices.push(baseIndex, baseIndex + 1, baseIndex + 2);
      indices.push(baseIndex + 3, baseIndex + 4, baseIndex + 5);
      indices.push(baseIndex, baseIndex + 1, baseIndex + 3, baseIndex + 1, baseIndex + 4, baseIndex + 3);
      indices.push(baseIndex + 1, baseIndex + 2, baseIndex + 4, baseIndex + 2, baseIndex + 5, baseIndex + 4);
      indices.push(baseIndex + 2, baseIndex, baseIndex + 5, baseIndex, baseIndex + 3, baseIndex + 5);
  
      baseIndex += 6;
    }
  
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.setIndex(indices);
    geometry.computeVertexNormals();
  
    const material = new THREE.MeshStandardMaterial({
      vertexColors: true,
      side: THREE.DoubleSide,
      flatShading: true,
      roughness: 1,
      metalness: 0.1,
    });
  
    const mesh = new THREE.Mesh(geometry, material);
  
    const wireframeMaterial = new THREE.MeshBasicMaterial({ color: 0x747678, wireframe: true });
    const wireframeMesh = new THREE.Mesh(geometry, wireframeMaterial);
    wireframeMesh.scale.multiplyScalar(1.001);
  
    const group = new THREE.Group();
    group.add(mesh);
    group.add(wireframeMesh);
  
    mesh.castShadow = true;
    mesh.receiveShadow = true;
  
    group.updateMatrix();
  
    return group;
  }
  


///////////////////////

// Helper function to create a side face (quad face between base and top)
function createSideFace(p1, p2, p1Extruded, p2Extruded) {
  const sideGeometry = new THREE.Geometry();
  sideGeometry.vertices.push(p1, p2, p2Extruded, p1Extruded);
  sideGeometry.faces.push(new THREE.Face3(0, 1, 2));  // Side face 1
  sideGeometry.faces.push(new THREE.Face3(0, 2, 3));  // Side face 2
  return sideGeometry;
}





export default App;