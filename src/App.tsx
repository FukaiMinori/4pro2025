//import './App.css'
import * as THREE from 'three'
import { Canvas, extend } from '@react-three/fiber'
import { Hud, OrthographicCamera, shaderMaterial } from '@react-three/drei'
import { useControls } from 'leva'
//import './styles.css'
import vertexShader from "./vshader.vert?raw";
//import fragmentShader from "./sankaku.frag?raw";
import fragmentShader from "./sphericalFibo9.frag?raw";
import { useEffect, useState } from 'react';
//import type { alphaT } from 'three/tsl';


declare global
{ namespace JSX
 { interface IntrinsicElements
        { "diffMaterial": any}
 }
}


function App() {
  const { radius, num, an, embedSmall, embedSmall2} = useControls({
    radius :{value:1.0, min: 0.5, max: 2.0, step:0.05},
    num :{value:130.0, min: 10.0, max: 1000.0, step:1.0},
    an :{value:0.0, min:-3.14159 / 2.0, max:3.14159 / 2.0, step:0.1},
    embedSmall :{value:1.7, min: 0.5, max: 2.0, step:0.005},
    embedSmall2 :{value:1.7, min: 0.5, max: 2.0, step:0.005}
  });

  const [uTime, setUTime] = useState(0.0);
  useEffect(()=>{
    const interval = setInterval(()=>{
      setUTime((uTime)=>uTime+0.1);
      //console.log(uTime);
    },100)
    return ()=>clearInterval(interval);
  },[uTime]
    )
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
			<diffMaterial key={DiffMaterial.key} glslVersion={THREE.GLSL3} radius={radius} num={num} an={an} uTime={uTime} embedSmall={embedSmall} embedSmall2={embedSmall2}/>
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
   num: 130.0,
   an: 0.0,
   embedSmall: 1.7,
   embedSmall2: 1.7,
   uTime: 0.0
  },
  vertexShader,
  fragmentShader
);
extend({DiffMaterial});

export default App
