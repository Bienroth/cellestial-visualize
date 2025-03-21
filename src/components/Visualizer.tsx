
import React, { useRef, useEffect, useState, useMemo } from 'react';
import * as THREE from 'three';
import { SpatialData, Point, ViewState, VisualizationMode } from '@/types/data';
import { getColorForValue, getColorForCluster } from '@/utils/colorMappings';
import ColorScale from './ColorScale';
import { geneExpressionColorScale, clusterColors } from '@/utils/colorMappings';

interface VisualizerProps {
  data: SpatialData | null;
  selectedGene?: string;
  mode: VisualizationMode;
  className?: string;
  viewState: ViewState;
  onViewStateChange: (viewState: ViewState) => void;
}

const Visualizer: React.FC<VisualizerProps> = ({
  data,
  selectedGene,
  mode,
  className,
  viewState,
  onViewStateChange
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const pointsRef = useRef<THREE.Points | null>(null);
  const animationFrameRef = useRef<number>(0);
  const isDraggingRef = useRef<boolean>(false);
  const isPanningRef = useRef<boolean>(false);
  const previousMousePosition = useRef<{ x: number; y: number }>({ x: 0, y: 0 });
  
  const [colorScale, setColorScale] = useState<typeof geneExpressionColorScale>(geneExpressionColorScale);
  
  // Calculate the visible expression data
  const expressionData = useMemo(() => {
    if (!data || !selectedGene || mode !== 'gene') return null;
    
    const values = data.expressionData[selectedGene] || [];
    if (values.length === 0) return null;
    
    const min = Math.min(...values);
    const max = Math.max(...values);
    
    setColorScale({
      ...geneExpressionColorScale,
      min,
      max
    });
    
    return values;
  }, [data, selectedGene, mode]);

  // Initialize scene, camera, and renderer
  useEffect(() => {
    if (!containerRef.current) return;
    
    // Scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xffffff);
    sceneRef.current = scene;
    
    // Camera
    const aspectRatio = containerRef.current.clientWidth / containerRef.current.clientHeight;
    const camera = new THREE.PerspectiveCamera(60, aspectRatio, 0.1, 1000);
    camera.position.z = 2;
    cameraRef.current = camera;
    
    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(containerRef.current.clientWidth, containerRef.current.clientHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    containerRef.current.appendChild(renderer.domElement);
    rendererRef.current = renderer;
    
    // Add ambient light
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.7);
    scene.add(ambientLight);
    
    // Add directional light
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
    directionalLight.position.set(1, 1, 1);
    scene.add(directionalLight);
    
    // Handle resize
    const handleResize = () => {
      if (!containerRef.current || !rendererRef.current || !cameraRef.current) return;
      
      const width = containerRef.current.clientWidth;
      const height = containerRef.current.clientHeight;
      
      rendererRef.current.setSize(width, height);
      cameraRef.current.aspect = width / height;
      cameraRef.current.updateProjectionMatrix();
    };
    
    window.addEventListener('resize', handleResize);
    
    // Animation loop
    const animate = () => {
      if (!sceneRef.current || !cameraRef.current || !rendererRef.current) return;
      
      rendererRef.current.render(sceneRef.current, cameraRef.current);
      animationFrameRef.current = requestAnimationFrame(animate);
    };
    
    animate();
    
    return () => {
      window.removeEventListener('resize', handleResize);
      cancelAnimationFrame(animationFrameRef.current);
      
      if (rendererRef.current && containerRef.current) {
        containerRef.current.removeChild(rendererRef.current.domElement);
      }
      
      if (sceneRef.current) {
        // Clean up all objects from the scene
        while (sceneRef.current.children.length > 0) {
          const object = sceneRef.current.children[0];
          sceneRef.current.remove(object);
        }
      }
      
      rendererRef.current?.dispose();
    };
  }, []);

  // Update camera and controls when viewState changes
  useEffect(() => {
    if (!cameraRef.current) return;
    
    cameraRef.current.position.z = 2 / viewState.zoom;
    cameraRef.current.position.x = viewState.target[0];
    cameraRef.current.position.y = viewState.target[1];
    
    cameraRef.current.rotation.x = viewState.rotationX;
    cameraRef.current.rotation.y = viewState.rotationY;
    
    cameraRef.current.updateProjectionMatrix();
  }, [viewState]);

  // Create or update points geometry based on data
  useEffect(() => {
    if (!data || !sceneRef.current) return;
    
    // Remove existing points if any
    if (pointsRef.current) {
      sceneRef.current.remove(pointsRef.current);
      pointsRef.current.geometry.dispose();
      (pointsRef.current.material as THREE.Material).dispose();
    }
    
    // Create geometry for points
    const geometry = new THREE.BufferGeometry();
    
    // Create positions array
    const positions = new Float32Array(data.points.length * 3);
    const colors = new Float32Array(data.points.length * 3);
    
    // Set positions and colors
    data.points.forEach((point, i) => {
      positions[i * 3] = point.x;
      positions[i * 3 + 1] = point.y;
      positions[i * 3 + 2] = point.z || 0;
      
      let color;
      if (mode === 'gene' && expressionData && selectedGene) {
        const value = expressionData[i] || 0;
        color = new THREE.Color(getColorForValue(value, colorScale));
      } else if (mode === 'cluster' && point.clusterId) {
        color = new THREE.Color(getColorForCluster(point.clusterId));
      } else {
        color = new THREE.Color(0xcccccc);
      }
      
      colors[i * 3] = color.r;
      colors[i * 3 + 1] = color.g;
      colors[i * 3 + 2] = color.b;
    });
    
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    
    // Create material with adjusted size for better performance with large datasets
    const material = new THREE.PointsMaterial({
      size: data.points.length > 100000 ? 0.005 : 0.01,
      vertexColors: true,
      sizeAttenuation: true,
      transparent: true,
      opacity: 0.8
    });
    
    // Create points
    const points = new THREE.Points(geometry, material);
    sceneRef.current.add(points);
    pointsRef.current = points;
    
  }, [data, mode, selectedGene, expressionData, colorScale]);

  // Mouse and touch event handlers for interaction
  useEffect(() => {
    if (!containerRef.current) return;
    
    const handleMouseDown = (e: MouseEvent) => {
      if (e.shiftKey) {
        // Panning mode when shift is pressed
        isPanningRef.current = true;
      } else {
        // Rotation mode by default
        isDraggingRef.current = true;
      }
      previousMousePosition.current = {
        x: e.clientX,
        y: e.clientY
      };
    };
    
    const handleMouseMove = (e: MouseEvent) => {
      if (isDraggingRef.current) {
        // Handle rotation
        const deltaMove = {
          x: e.clientX - previousMousePosition.current.x,
          y: e.clientY - previousMousePosition.current.y
        };
        
        onViewStateChange({
          ...viewState,
          rotationY: viewState.rotationY + deltaMove.x * 0.01,
          rotationX: viewState.rotationX + deltaMove.y * 0.01
        });
      } else if (isPanningRef.current) {
        // Handle panning (left, right, front, back)
        const deltaMove = {
          x: e.clientX - previousMousePosition.current.x,
          y: e.clientY - previousMousePosition.current.y
        };
        
        // Scale the movement based on zoom level
        const panFactor = 0.005 / viewState.zoom;
        
        onViewStateChange({
          ...viewState,
          target: [
            viewState.target[0] - deltaMove.x * panFactor,
            viewState.target[1] + deltaMove.y * panFactor,
            viewState.target[2]
          ]
        });
      }
      
      previousMousePosition.current = {
        x: e.clientX,
        y: e.clientY
      };
    };
    
    const handleMouseUp = () => {
      isDraggingRef.current = false;
      isPanningRef.current = false;
    };
    
    const handleWheel = (e: WheelEvent) => {
      e.preventDefault();
      
      const zoomDelta = e.deltaY * -0.001;
      const newZoom = Math.max(0.5, Math.min(5, viewState.zoom + zoomDelta));
      
      onViewStateChange({
        ...viewState,
        zoom: newZoom
      });
    };
    
    const handleKeyDown = (e: KeyboardEvent) => {
      // Pan with arrow keys
      const panAmount = 0.05 / viewState.zoom;
      let newTarget = [...viewState.target];
      
      switch (e.key) {
        case 'ArrowLeft':
          newTarget[0] += panAmount;
          break;
        case 'ArrowRight':
          newTarget[0] -= panAmount;
          break;
        case 'ArrowUp':
          newTarget[1] -= panAmount;
          break;
        case 'ArrowDown':
          newTarget[1] += panAmount;
          break;
      }
      
      if (e.key.startsWith('Arrow')) {
        onViewStateChange({
          ...viewState,
          target: newTarget as [number, number, number]
        });
      }
    };
    
    const container = containerRef.current;
    container.addEventListener('mousedown', handleMouseDown);
    window.addEventListener('mousemove', handleMouseMove);
    window.addEventListener('mouseup', handleMouseUp);
    container.addEventListener('wheel', handleWheel, { passive: false });
    window.addEventListener('keydown', handleKeyDown);
    
    return () => {
      container.removeEventListener('mousedown', handleMouseDown);
      window.removeEventListener('mousemove', handleMouseMove);
      window.removeEventListener('mouseup', handleMouseUp);
      container.removeEventListener('wheel', handleWheel);
      window.removeEventListener('keydown', handleKeyDown);
    };
  }, [viewState, onViewStateChange]);

  // Generate color legend
  const renderLegend = () => {
    if (mode === 'gene' && selectedGene) {
      return (
        <div className="absolute bottom-4 right-4 w-64 glass p-3 rounded-lg">
          <ColorScale
            scale={colorScale}
            label={`${selectedGene} Expression`}
          />
        </div>
      );
    }
    
    if (mode === 'cluster' && data?.clusters) {
      return (
        <div className="absolute bottom-4 right-4 max-w-sm max-h-60 glass p-3 rounded-lg overflow-auto">
          <div className="text-xs font-medium text-muted-foreground mb-2">Clusters</div>
          <div className="grid grid-cols-2 gap-1">
            {Object.entries(data.clusters).map(([id, name]) => (
              <div key={id} className="flex items-center text-xs space-x-1.5">
                <div
                  className="w-3 h-3 rounded-full flex-shrink-0"
                  style={{ backgroundColor: getColorForCluster(id) }}
                />
                <span className="truncate">{name}</span>
              </div>
            ))}
          </div>
        </div>
      );
    }
    
    return null;
  };

  return (
    <div className={`relative ${className}`}>
      <div ref={containerRef} className="w-full h-full" />
      {renderLegend()}
      
      {/* Controls help tooltip */}
      <div className="absolute bottom-4 left-4 glass p-2 rounded-lg text-xs">
        <p className="font-medium">Controls:</p>
        <ul className="text-muted-foreground space-y-1 mt-1">
          <li>Drag: Rotate view</li>
          <li>Shift+Drag: Pan view</li>
          <li>Arrow keys: Pan view</li>
          <li>Scroll: Zoom in/out</li>
        </ul>
      </div>
    </div>
  );
};

export default Visualizer;
