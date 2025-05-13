//import './App.css'
import * as THREE from 'three'
import { Canvas, extend } from '@react-three/fiber'
import { Hud, OrthographicCamera, shaderMaterial } from '@react-three/drei'
import { useControls } from 'leva'
//import './styles.css'
import vertexShader from "./vshader.vert?raw";
//import fragmentShader from "./sankaku.frag?raw";
import fragmentShader from "./may13.frag?raw";
import type { alphaT } from 'three/tsl';


declare global
{ namespace JSX
 { interface IntrinsicElements
        { "diffMaterial": any}
 }
}


function App() {
  const { radius, alpha } = useControls({radius :{value:1.0, min: 0.5, max: 2.0, step:0.05},
					alpha :{value:0.7, min: 0.0, max: 1.0, step:0.05}
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
			<diffMaterial key={DiffMaterial.key} glslVersion={THREE.GLSL3} radius={radius} alpha={alpha}/>
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
   alpha: 1.0
  },
  vertexShader,
  fragmentShader
);
extend({DiffMaterial});

export default App
