import React from 'react';

const PlaceholderMolecule = () => (
  <svg 
    width="100%" 
    height="100%" 
    viewBox="0 0 200 200" 
    xmlns="http://www.w3.org/2000/svg"
  >
    <circle cx="100" cy="100" r="50" fill="#f0f0f0" stroke="#aaaaaa" strokeWidth="2" />
    <circle cx="100" cy="70" r="12" fill="#aaaaaa" />
    <circle cx="70" cy="115" r="12" fill="#aaaaaa" />
    <circle cx="130" cy="115" r="12" fill="#aaaaaa" />
    <line x1="100" y1="70" x2="70" y2="115" stroke="#aaaaaa" strokeWidth="3" />
    <line x1="100" y1="70" x2="130" y2="115" stroke="#aaaaaa" strokeWidth="3" />
    <line x1="70" y1="115" x2="130" y2="115" stroke="#aaaaaa" strokeWidth="3" />
  </svg>
);

export default PlaceholderMolecule;