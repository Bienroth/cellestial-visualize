
import React from 'react';
import { Button } from '@/components/ui/button';
import { Slider } from '@/components/ui/slider';
import {
  ZoomIn,
  ZoomOut,
  RotateCcw,
  Home,
  Download,
  Maximize2,
  Minimize2
} from 'lucide-react';
import { ViewState } from '@/types/data';

interface ControlPanelProps {
  viewState: ViewState;
  onZoomIn: () => void;
  onZoomOut: () => void;
  onReset: () => void;
  onDownload: () => void;
  onFullscreen: () => void;
  isFullscreen: boolean;
  className?: string;
}

const ControlPanel: React.FC<ControlPanelProps> = ({
  viewState,
  onZoomIn,
  onZoomOut,
  onReset,
  onDownload,
  onFullscreen,
  isFullscreen,
  className
}) => {
  return (
    <div className={`glass rounded-lg shadow-sm backdrop-blur-md ${className}`}>
      <div className="p-2 grid grid-cols-1 gap-1">
        <Button
          variant="ghost"
          size="icon"
          onClick={onZoomIn}
          className="h-8 w-8"
          title="Zoom In"
        >
          <ZoomIn className="h-4 w-4" />
        </Button>
        
        <Button
          variant="ghost"
          size="icon"
          onClick={onZoomOut}
          className="h-8 w-8"
          title="Zoom Out"
        >
          <ZoomOut className="h-4 w-4" />
        </Button>
        
        <div className="my-1 w-full h-px bg-border" />
        
        <Button
          variant="ghost"
          size="icon"
          onClick={onReset}
          className="h-8 w-8"
          title="Reset View"
        >
          <Home className="h-4 w-4" />
        </Button>
        
        <Button
          variant="ghost"
          size="icon"
          onClick={onDownload}
          className="h-8 w-8"
          title="Download Image"
        >
          <Download className="h-4 w-4" />
        </Button>
        
        <Button
          variant="ghost"
          size="icon"
          onClick={onFullscreen}
          className="h-8 w-8"
          title={isFullscreen ? "Exit Fullscreen" : "Enter Fullscreen"}
        >
          {isFullscreen ? 
            <Minimize2 className="h-4 w-4" /> : 
            <Maximize2 className="h-4 w-4" />
          }
        </Button>
      </div>
    </div>
  );
};

export default ControlPanel;
