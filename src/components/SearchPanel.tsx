
import React, { useState, useEffect, useRef } from 'react';
import { Input } from '@/components/ui/input';
import { Button } from '@/components/ui/button';
import { ScrollArea } from '@/components/ui/scroll-area';
import { 
  Search,
  X,
  ArrowRight,
  Layers
} from 'lucide-react';

interface SearchPanelProps {
  genes: string[];
  onGeneSelect: (gene: string) => void;
  onToggleClusterView: () => void;
  isClusterViewActive: boolean;
  selectedGene?: string;
  className?: string;
}

const SearchPanel: React.FC<SearchPanelProps> = ({
  genes,
  onGeneSelect,
  onToggleClusterView,
  isClusterViewActive,
  selectedGene,
  className,
}) => {
  const [searchQuery, setSearchQuery] = useState('');
  const [filteredGenes, setFilteredGenes] = useState<string[]>([]);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (searchQuery.trim() === '') {
      setFilteredGenes([]);
    } else {
      const query = searchQuery.toLowerCase();
      const filtered = genes
        .filter(gene => gene.toLowerCase().includes(query))
        .slice(0, 50); // Limit results for performance
      setFilteredGenes(filtered);
    }
  }, [searchQuery, genes]);

  const handleClear = () => {
    setSearchQuery('');
    inputRef.current?.focus();
  };

  return (
    <div className={`flex flex-col h-full bg-card border-l border-border ${className}`}>
      <div className="p-4 border-b">
        <h3 className="text-lg font-semibold mb-1">CELLestial</h3>
        <p className="text-sm text-muted-foreground">
          Spatial Transcriptomics Visualizer
        </p>
      </div>
      
      <div className="p-4 border-b">
        <div className="space-y-4">
          <div>
            <h4 className="text-sm font-medium mb-1.5">Gene Expression</h4>
            <div className="relative">
              <div className="absolute inset-y-0 left-0 flex items-center pl-2.5 pointer-events-none">
                <Search className="h-4 w-4 text-muted-foreground" />
              </div>
              <Input
                ref={inputRef}
                type="text"
                placeholder="Search genes..."
                value={searchQuery}
                onChange={e => setSearchQuery(e.target.value)}
                className="pl-8 pr-8"
              />
              {searchQuery && (
                <button
                  onClick={handleClear}
                  className="absolute inset-y-0 right-0 flex items-center pr-2.5"
                >
                  <X className="h-4 w-4 text-muted-foreground hover:text-foreground" />
                </button>
              )}
            </div>
          </div>
          
          <Button 
            variant={isClusterViewActive ? "default" : "outline"}
            className="w-full flex items-center justify-between"
            onClick={onToggleClusterView}
          >
            <span className="flex items-center">
              <Layers className="mr-2 h-4 w-4" />
              Visualize Clusters
            </span>
            <ArrowRight className="h-4 w-4" />
          </Button>
        </div>
      </div>
      
      {filteredGenes.length > 0 && (
        <div className="flex-1 overflow-hidden border-b">
          <ScrollArea className="h-full">
            <div className="p-2">
              <div className="text-xs font-medium text-muted-foreground mb-1 px-2">
                {filteredGenes.length} genes found
              </div>
              <div className="space-y-1">
                {filteredGenes.map(gene => (
                  <button
                    key={gene}
                    onClick={() => {
                      onGeneSelect(gene);
                      setSearchQuery('');
                    }}
                    className={`w-full text-left px-3 py-1.5 rounded-md text-sm transition-colors
                      ${selectedGene === gene 
                        ? 'bg-primary text-primary-foreground'
                        : 'hover:bg-accent hover:text-accent-foreground'
                      }`}
                  >
                    {gene}
                  </button>
                ))}
              </div>
            </div>
          </ScrollArea>
        </div>
      )}
      
      {selectedGene && !isClusterViewActive && (
        <div className="p-4">
          <div className="text-sm font-medium mb-1">Selected Gene</div>
          <div className="flex items-center justify-between bg-accent rounded-md p-2">
            <span className="text-sm">{selectedGene}</span>
            <Button 
              variant="ghost" 
              size="sm"
              onClick={() => onGeneSelect('')}
              className="h-6 w-6 p-0"
            >
              <X className="h-3 w-3" />
            </Button>
          </div>
        </div>
      )}
      
      <div className="p-4 mt-auto">
        <div className="text-xs text-muted-foreground">
          Stereo-seq spatial transcriptomics visualization
        </div>
      </div>
    </div>
  );
};

export default SearchPanel;
