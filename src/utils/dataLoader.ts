
import { SpatialData, Point } from '@/types/data';

/**
 * Mock function to simulate loading h5ad data
 * In a real implementation, this would use a library like hdf5.js or a backend service
 */
export const loadH5adData = async (file: File): Promise<SpatialData> => {
  // This is a mock implementation
  return new Promise((resolve) => {
    console.log(`Processing file: ${file.name}`);
    
    // Simulate loading delay
    setTimeout(() => {
      // For demonstration, we'll create mock data
      // In a real implementation, this would parse the h5ad file
      const mockData = createMockData();
      resolve(mockData);
    }, 1500);
  });
};

/**
 * Creates mock data for development and testing
 */
function createMockData(): SpatialData {
  // Create a grid of points
  const points: Point[] = [];
  const size = 100;
  const geneCount = 500;
  
  // Generate points in a grid pattern with some noise
  for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
      if (Math.random() > 0.7) { // Add some sparsity
        const noise = Math.random() * 0.3;
        points.push({
          x: i / size * 2 - 1 + noise,
          y: j / size * 2 - 1 + noise,
          z: Math.random() * 0.2,
          clusterId: String(Math.floor(Math.random() * 10))
        });
      }
    }
  }
  
  // Generate gene names
  const genes = Array.from({ length: geneCount }, (_, i) => `Gene${i + 1}`);
  
  // Generate mock expression data
  const expressionData: { [gene: string]: number[] } = {};
  genes.forEach(gene => {
    expressionData[gene] = points.map(() => Math.random() * 10);
  });
  
  // Generate cluster names
  const clusters: { [key: string]: string } = {};
  for (let i = 0; i < 10; i++) {
    clusters[i.toString()] = `Cluster ${i + 1}`;
  }
  
  return {
    points,
    genes,
    clusters,
    expressionData
  };
}

/**
 * Normalizes expression values to a range of [0, 1]
 */
export const normalizeExpressionValues = (
  values: number[]
): { normalizedValues: number[], min: number, max: number } => {
  const min = Math.min(...values);
  const max = Math.max(...values);
  const range = max - min;
  
  const normalizedValues = values.map(val => {
    if (range === 0) return 0;
    return (val - min) / range;
  });
  
  return { normalizedValues, min, max };
};
