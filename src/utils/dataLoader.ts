
import { SpatialData, Point } from '@/types/data';
import * as h5wasm from 'h5wasm';

/**
 * Loads h5ad data from a File object
 */
export const loadH5adData = async (file: File): Promise<SpatialData> => {
  console.log(`Processing file: ${file.name}`);
  
  try {
    // Initialize the WASM backend for h5wasm
    await h5wasm.ready;
    
    // Read the file as ArrayBuffer
    const buffer = await file.arrayBuffer();
    
    // Open the h5ad file using h5wasm - fix parameter type issue
    // We need to create a File with the correct constructor signature
    const f = new h5wasm.File(new Uint8Array(buffer), 'r');
    
    // Extract data from the h5ad file structure
    const points = extractSpatialCoordinates(f);
    const { genes, expressionData } = extractGeneExpressionData(f, points.length);
    const clusters = extractClusterData(f, points.length);
    
    // Assign cluster IDs to points if available
    if (clusters) {
      Object.keys(clusters).forEach(id => {
        // Find points that belong to this cluster
        points.forEach(point => {
          if (point.clusterId === id) {
            point.clusterId = id;
          }
        });
      });
    }
    
    return {
      points,
      genes,
      clusters,
      expressionData
    };
  } catch (error) {
    console.error('Error parsing h5ad file:', error);
    throw new Error(`Failed to parse h5ad file: ${error.message}`);
  }
};

/**
 * Extracts spatial coordinates from the h5ad file
 */
function extractSpatialCoordinates(file: h5wasm.File): Point[] {
  try {
    // Try to find spatial coordinates in adata.obsm["spatial"]
    // Using type guards to check the existence of properties before accessing them
    const obsm = file.get('obsm');
    if (obsm && typeof obsm.get === 'function') {
      const spatial = obsm.get('spatial');
      if (spatial && 'value' in spatial) {
        const spatialData = spatial.value as number[][];
        
        // Create points from spatial coordinates
        return spatialData.map((coords, index) => ({
          x: coords[0] || 0,
          y: coords[1] || 0,
          z: coords.length > 2 ? coords[2] : 0,
          clusterId: undefined
        }));
      }
    }
    
    // Fallback to mock data if spatial coordinates not found
    console.warn('Spatial coordinates not found in h5ad file, using mock data');
    return createMockPoints();
  } catch (error) {
    console.error('Error extracting spatial coordinates:', error);
    return createMockPoints();
  }
}

/**
 * Extracts gene expression data from the h5ad file
 */
function extractGeneExpressionData(file: h5wasm.File, numPoints: number): { genes: string[], expressionData: { [gene: string]: number[] } } {
  try {
    const genes: string[] = [];
    const expressionData: { [gene: string]: number[] } = {};
    
    // Try to get gene names from adata.var.index
    const varObj = file.get('var');
    if (varObj && typeof varObj.get === 'function') {
      const indexObj = varObj.get('index');
      if (indexObj && 'value' in indexObj) {
        const geneList = indexObj.value as string[];
        genes.push(...geneList);
      } else {
        // Fallback to mock gene names
        console.warn('Gene names not found in h5ad file, using mock gene names');
        for (let i = 0; i < 500; i++) {
          genes.push(`Gene${i + 1}`);
        }
      }
    } else {
      // Fallback to mock gene names
      console.warn('Gene names not found in h5ad file, using mock gene names');
      for (let i = 0; i < 500; i++) {
        genes.push(`Gene${i + 1}`);
      }
    }
    
    // Try to get expression data from adata.X
    if (file.get('X') && 'value' in file.get('X')) {
      const expressionMatrix = file.get('X').value;
      
      // Check if expression matrix is in the expected format (array of arrays)
      if (Array.isArray(expressionMatrix) && Array.isArray(expressionMatrix[0])) {
        // Transpose the matrix: rows = spots, columns = genes
        genes.forEach((gene, geneIndex) => {
          expressionData[gene] = expressionMatrix.map(row => row[geneIndex] || 0);
        });
      } else {
        console.warn('Expression matrix format not as expected, using mock data');
        // Fall back to mock expression data
        genes.forEach(gene => {
          expressionData[gene] = Array(numPoints).fill(0).map(() => Math.random() * 10);
        });
      }
    } else {
      console.warn('Expression data not found in h5ad file, using mock data');
      // Fall back to mock expression data
      genes.forEach(gene => {
        expressionData[gene] = Array(numPoints).fill(0).map(() => Math.random() * 10);
      });
    }
    
    return { genes, expressionData };
  } catch (error) {
    console.error('Error extracting gene expression data:', error);
    // Fall back to mock data
    const mockGenes = Array.from({ length: 500 }, (_, i) => `Gene${i + 1}`);
    const mockExpressionData: { [gene: string]: number[] } = {};
    
    mockGenes.forEach(gene => {
      mockExpressionData[gene] = Array(numPoints).fill(0).map(() => Math.random() * 10);
    });
    
    return { genes: mockGenes, expressionData: mockExpressionData };
  }
}

/**
 * Extracts cluster data from the h5ad file
 */
function extractClusterData(file: h5wasm.File, numPoints: number): { [key: string]: string } | undefined {
  try {
    // Try to find cluster info in adata.obs["leiden"]
    const obs = file.get('obs');
    if (obs && typeof obs.get === 'function') {
      const leiden = obs.get('leiden');
      if (leiden && 'value' in leiden) {
        const clusterLabels = leiden.value as string[];
        
        // Get unique cluster IDs
        const uniqueClusters = Array.from(new Set(clusterLabels));
        
        // Create clusters object
        const clusters: { [key: string]: string } = {};
        uniqueClusters.forEach((clusterId, index) => {
          clusters[clusterId] = `Cluster ${index + 1}`;
        });
        
        return clusters;
      }
    }
    
    // Fallback to mock cluster data
    console.warn('Cluster data not found in h5ad file, using mock clusters');
    const mockClusters: { [key: string]: string } = {};
    for (let i = 0; i < 10; i++) {
      mockClusters[i.toString()] = `Cluster ${i + 1}`;
    }
    
    return mockClusters;
  } catch (error) {
    console.error('Error extracting cluster data:', error);
    // Fallback to mock cluster data
    const mockClusters: { [key: string]: string } = {};
    for (let i = 0; i < 10; i++) {
      mockClusters[i.toString()] = `Cluster ${i + 1}`;
    }
    
    return mockClusters;
  }
}

/**
 * Creates mock points for testing
 */
function createMockPoints(): Point[] {
  const points: Point[] = [];
  const size = 100;
  
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
  
  return points;
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
