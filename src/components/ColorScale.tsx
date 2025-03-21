
import React from 'react';
import { ColorScale as ColorScaleType } from '@/types/data';
import { generateGradientSteps } from '@/utils/colorMappings';

interface ColorScaleProps {
  scale: ColorScaleType;
  label: string;
  className?: string;
}

const ColorScale: React.FC<ColorScaleProps> = ({ scale, label, className }) => {
  const gradientSteps = generateGradientSteps(scale, 10);
  
  return (
    <div className={`flex flex-col space-y-1 ${className}`}>
      <div className="text-xs font-medium text-muted-foreground">{label}</div>
      <div className="flex items-center space-x-1">
        <div className="w-full h-3 rounded-md overflow-hidden flex">
          {gradientSteps.map((color, i) => (
            <div 
              key={i} 
              className="h-full flex-1" 
              style={{ backgroundColor: color }}
            />
          ))}
        </div>
      </div>
      <div className="flex justify-between text-xs text-muted-foreground">
        <span>{scale.min.toFixed(1)}</span>
        <span>{scale.max.toFixed(1)}</span>
      </div>
    </div>
  );
};

export default ColorScale;
