
import { ColorScale } from '@/types/data';

// Viridis-inspired color palette for gene expression
export const geneExpressionColorScale: ColorScale = {
  min: 0,
  max: 10,
  colors: [
    '#440154', '#472c7a', '#3b518b', '#2c718e', '#21908d', 
    '#27ad81', '#5cc863', '#aadc32', '#fde725'
  ]
};

// Color palette for clusters (distinct colors)
export const clusterColors = [
  '#4e79a7', '#f28e2c', '#e15759', '#76b7b2', '#59a14f',
  '#edc949', '#af7aa1', '#ff9da7', '#9c755f', '#bab0ab'
];

/**
 * Get color for a value based on a color scale
 */
export const getColorForValue = (
  value: number, 
  scale: ColorScale = geneExpressionColorScale
): string => {
  // Handle edge cases
  if (isNaN(value) || value === undefined) return '#cccccc';
  
  // Normalize value to [0, 1] range based on the scale
  const normalizedValue = Math.max(0, Math.min(1, (value - scale.min) / (scale.max - scale.min)));
  
  // Find color based on the normalized value
  const colorIndex = Math.min(
    scale.colors.length - 1,
    Math.floor(normalizedValue * scale.colors.length)
  );
  
  return scale.colors[colorIndex];
};

/**
 * Get color for a cluster ID
 */
export const getColorForCluster = (clusterId: string): string => {
  const numericId = parseInt(clusterId, 10);
  if (isNaN(numericId)) return '#cccccc';
  
  return clusterColors[numericId % clusterColors.length];
};

/**
 * Generate a gradient array for the color legend
 */
export const generateGradientSteps = (scale: ColorScale, steps: number = 10): string[] => {
  const result: string[] = [];
  for (let i = 0; i < steps; i++) {
    const value = scale.min + (i / (steps - 1)) * (scale.max - scale.min);
    result.push(getColorForValue(value, scale));
  }
  return result;
};
