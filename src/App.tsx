//import './App.css'
import * as THREE from 'three'
import { Canvas, extend } from '@react-three/fiber'
import { Hud, OrthographicCamera, shaderMaterial } from '@react-three/drei'
import { useControls } from 'leva'
//import './styles.css'
import vertexShader from "./vshader.vert?raw";
//import fragmentShader from "./sankaku.frag?raw";
import fragmentShader from "./may13.frag?raw";
//import type { alphaT } from 'three/tsl';


declare global
{ namespace JSX
 { interface IntrinsicElements
        { "diffMaterial": any}
 }
}


function App() {
  const { radius, wallZ, theta } = useControls({radius :{value:1.0, min: 0.5, max: 2.0, step:0.05},
					wallZ :{value:0.0, min: -0.15, max: 0.3, step:0.005},
          theta :{value:0.0, min: -3.14159*0.1, max: 3.14159*0.1, step:0.005}
				       });
  return (
    <>
	<div style={{ height: "100dvh", width: "100dvw" }}>
	    <Canvas>
		<Hud>
		    <OrthographicCamera
			makeDefault
			top={1}
			right={1}
			bottom={-1}
			left={-1}
			near={0}
			far={1}
		    />
		    <OrthographicCamera  makeDefault top={1} right={1} bottom={-1} left={-1} near={0} far={1} />
		    <mesh>
			<planeGeometry args={[2,2]}/>
			{/* @ts-ignore TS2339: Property 'diffMaterial' does not exist on type 'JSX.IntrinsicElements'.*/}
			<diffMaterial key={DiffMaterial.key} glslVersion={THREE.GLSL3} radius={radius} wallZ={wallZ} theta={theta}/>
		    </mesh>
		</Hud>
	    </Canvas>
	</div>
    </>
  )
}

const DiffMaterial = shaderMaterial(
  {u_resolution: new THREE.Vector2(window.innerWidth, window.innerHeight),
   radius: 1.0,
   wallZ: 0.0,
   theta: 0.0
  },
  vertexShader,
  fragmentShader
);
extend({DiffMaterial});

export default App
