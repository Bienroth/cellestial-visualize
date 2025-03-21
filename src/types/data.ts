
export interface SpatialData {
  points: Point[];
  genes: string[];
  clusters?: { [key: string]: string };
  expressionData: { [gene: string]: number[] };
  // Additional metadata can be added here as we extract more info from h5ad
}

export interface Point {
  x: number;
  y: number;
  z?: number;
  clusterId?: string;
}

export interface ViewState {
  zoom: number;
  target: [number, number, number];
  rotationX: number;
  rotationY: number;
}

export interface ColorScale {
  min: number;
  max: number;
  colors: string[];
}

export type VisualizationMode = 'gene' | 'cluster';

export interface SelectedData {
  mode: VisualizationMode;
  geneId?: string;
  expressionRange?: [number, number];
}

export interface H5adMetadata {
  sampleId?: string;
  date?: string;
  experimentType?: string;
  description?: string;
}
