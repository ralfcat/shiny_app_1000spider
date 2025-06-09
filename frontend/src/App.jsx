import React, { useState } from 'react';

function App() {
  const [sequence, setSequence] = useState('');
  const [predictions, setPredictions] = useState(null);
  const [loading, setLoading] = useState(false);

  const submitSequence = async () => {
    setLoading(true);
    try {
      const res = await fetch('/predict', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence })
      });
      const data = await res.json();
      setPredictions(data);
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container">
      <h1>Spider Gene Visualizer</h1>
      <textarea value={sequence} onChange={(e) => setSequence(e.target.value)} placeholder="Enter FASTA sequence" rows="5" />
      <button onClick={submitSequence} disabled={loading}>
        {loading ? 'Loading...' : 'Predict Expression'}
      </button>
      {predictions && (
        <pre>{JSON.stringify(predictions, null, 2)}</pre>
      )}
    </div>
  );
}

export default App;
