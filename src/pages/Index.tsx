
import React, { useState, useEffect, useCallback } from 'react';
import Visualizer from '@/components/Visualizer';
import SearchPanel from '@/components/SearchPanel';
import ControlPanel from '@/components/ControlPanel';
import { Button } from '@/components/ui/button';
import { useToast } from '@/components/ui/use-toast';
import { SpatialData, ViewState, VisualizationMode } from '@/types/data';
import { loadH5adData, createMockData } from '@/utils/dataLoader';
import { Upload, FileUp, Play } from 'lucide-react';

const Index = () => {
  const { toast } = useToast();
  const [loading, setLoading] = useState(false);
  const [data, setData] = useState<SpatialData | null>(null);
  const [selectedGene, setSelectedGene] = useState<string>('');
  const [mode, setMode] = useState<VisualizationMode>('gene');
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [isPanelOpen, setIsPanelOpen] = useState(true);
  
  const [viewState, setViewState] = useState<ViewState>({
    zoom: 1,
    target: [0, 0, 0],
    rotationX: 0,
    rotationY: 0
  });

  // Handle file upload
  const handleFileUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;
    
    try {
      setLoading(true);
      
      const loadedData = await loadH5adData(file);
      setData(loadedData);
      
      toast({
        title: "Data loaded successfully",
        description: `Loaded ${loadedData.points.length} data points from ${file.name}`,
      });
      
    } catch (error) {
      console.error('Error loading data:', error);
      toast({
        title: "Error loading data",
        description: "There was an error loading the h5ad file.",
        variant: "destructive",
      });
    } finally {
      setLoading(false);
    }
  };

  // Handle loading demo data
  const handleLoadDemoData = () => {
    try {
      setLoading(true);
      
      // Use the mock data generator
      const mockData = createMockData();
      setData(mockData);
      
      toast({
        title: "Demo data loaded",
        description: `Loaded ${mockData.points.length} sample data points`,
      });
    } catch (error) {
      console.error('Error loading demo data:', error);
      toast({
        title: "Error loading demo data",
        description: "There was an error creating the demo data.",
        variant: "destructive",
      });
    } finally {
      setLoading(false);
    }
  };

  // Handle gene selection
  const handleGeneSelect = (gene: string) => {
    setSelectedGene(gene);
    if (gene) {
      setMode('gene');
    }
  };

  // Toggle between gene expression and cluster view
  const handleToggleClusterView = () => {
    setMode(prev => prev === 'gene' ? 'cluster' : 'gene');
    if (mode === 'gene') {
      // When switching to cluster view, clear gene selection
      setSelectedGene('');
    }
  };

  // Control panel actions
  const handleZoomIn = useCallback(() => {
    setViewState(prev => ({
      ...prev,
      zoom: Math.min(5, prev.zoom * 1.2)
    }));
  }, []);

  const handleZoomOut = useCallback(() => {
    setViewState(prev => ({
      ...prev,
      zoom: Math.max(0.5, prev.zoom / 1.2)
    }));
  }, []);

  const handleResetView = useCallback(() => {
    setViewState({
      zoom: 1,
      target: [0, 0, 0],
      rotationX: 0,
      rotationY: 0
    });
  }, []);

  // Download current view as image
  const handleDownload = useCallback(() => {
    const canvas = document.querySelector('canvas');
    if (!canvas) return;
    
    // Create a download link
    const link = document.createElement('a');
    link.download = `cellestial-${mode === 'gene' ? selectedGene || 'viz' : 'clusters'}.png`;
    link.href = canvas.toDataURL('image/png');
    link.click();
    
    toast({
      title: "Image downloaded",
      description: `The current view has been saved as ${link.download}`,
    });
  }, [mode, selectedGene, toast]);

  // Toggle fullscreen
  const handleFullscreen = useCallback(() => {
    if (!document.fullscreenElement) {
      document.documentElement.requestFullscreen().catch(err => {
        console.error(`Error attempting to enable full-screen mode: ${err.message}`);
      });
      setIsFullscreen(true);
    } else {
      if (document.exitFullscreen) {
        document.exitFullscreen();
        setIsFullscreen(false);
      }
    }
  }, []);

  // Handle fullscreen change events
  useEffect(() => {
    const handleFullscreenChange = () => {
      setIsFullscreen(!!document.fullscreenElement);
    };
    
    document.addEventListener('fullscreenchange', handleFullscreenChange);
    return () => {
      document.removeEventListener('fullscreenchange', handleFullscreenChange);
    };
  }, []);

  return (
    <div className="flex flex-col h-screen bg-background">
      {/* Header */}
      <header className="border-b border-border bg-card py-2 px-4 flex items-center justify-between">
        <div className="flex items-center space-x-4">
          <h1 className="text-lg font-medium">CELLestial Visualizer</h1>
          {!data && (
            <div className="flex space-x-2">
              <label className="relative">
                <Button 
                  variant="outline" 
                  className="space-x-2 flex items-center"
                  disabled={loading}
                >
                  {loading ? (
                    <div className="spinner h-4 w-4" />
                  ) : (
                    <FileUp className="h-4 w-4" />
                  )}
                  <span>Load H5AD File</span>
                </Button>
                <input 
                  type="file" 
                  accept=".h5ad" 
                  className="absolute inset-0 w-full h-full opacity-0 cursor-pointer" 
                  onChange={handleFileUpload}
                  disabled={loading}
                />
              </label>
              
              <Button
                variant="secondary"
                className="space-x-2 flex items-center"
                onClick={handleLoadDemoData}
                disabled={loading}
              >
                {loading ? (
                  <div className="spinner h-4 w-4" />
                ) : (
                  <Play className="h-4 w-4" />
                )}
                <span>Demo</span>
              </Button>
            </div>
          )}
        </div>
        
        {data && (
          <div className="text-sm text-muted-foreground">
            {data.points.length.toLocaleString()} data points
          </div>
        )}
      </header>
      
      {/* Main content */}
      <main className="flex-1 flex overflow-hidden">
        {/* Visualization area */}
        <div className="relative flex-1 overflow-hidden bg-[#fafafa] dark:bg-background">
          {!data ? (
            <div className="h-full flex flex-col items-center justify-center">
              <div className="glass rounded-lg p-8 max-w-md text-center space-y-4 animate-scale-in">
                <h2 className="text-2xl font-semibold">Welcome to CELLestial</h2>
                <p className="text-muted-foreground">
                  Upload an h5ad file containing Stereo-seq data to start exploring spatial transcriptomics.
                </p>
                <div className="mt-4 flex flex-col sm:flex-row gap-2">
                  <label className="relative flex-grow">
                    <Button 
                      variant="default" 
                      size="lg"
                      className="w-full space-x-2"
                      disabled={loading}
                    >
                      {loading ? (
                        <div className="spinner h-5 w-5" />
                      ) : (
                        <Upload className="h-5 w-5" />
                      )}
                      <span>Upload H5AD File</span>
                    </Button>
                    <input 
                      type="file" 
                      accept=".h5ad" 
                      className="absolute inset-0 w-full h-full opacity-0 cursor-pointer" 
                      onChange={handleFileUpload}
                      disabled={loading}
                    />
                  </label>
                  
                  <Button
                    variant="secondary"
                    size="lg"
                    className="w-full sm:w-auto space-x-2"
                    onClick={handleLoadDemoData}
                    disabled={loading}
                  >
                    {loading ? (
                      <div className="spinner h-5 w-5" />
                    ) : (
                      <Play className="h-5 w-5" />
                    )}
                    <span>View Demo</span>
                  </Button>
                </div>
              </div>
            </div>
          ) : (
            <>
              <Visualizer 
                data={data}
                selectedGene={selectedGene}
                mode={mode}
                viewState={viewState}
                onViewStateChange={setViewState}
                className="h-full w-full"
              />
              
              <div className="absolute top-4 right-4">
                <ControlPanel 
                  viewState={viewState}
                  onZoomIn={handleZoomIn}
                  onZoomOut={handleZoomOut}
                  onReset={handleResetView}
                  onDownload={handleDownload}
                  onFullscreen={handleFullscreen}
                  isFullscreen={isFullscreen}
                />
              </div>
              
              <Button
                variant="outline"
                size="sm"
                className="absolute top-4 left-4 glass"
                onClick={() => setIsPanelOpen(!isPanelOpen)}
              >
                {isPanelOpen ? "Hide Panel" : "Show Panel"}
              </Button>
            </>
          )}
        </div>
        
        {/* Side panel */}
        {data && isPanelOpen && (
          <div className="w-64 h-full flex-shrink-0 transition-all duration-300 ease-out-expo">
            <SearchPanel 
              genes={data.genes}
              onGeneSelect={handleGeneSelect}
              onToggleClusterView={handleToggleClusterView}
              isClusterViewActive={mode === 'cluster'}
              selectedGene={selectedGene}
              className="h-full"
            />
          </div>
        )}
      </main>
    </div>
  );
};

export default Index;
