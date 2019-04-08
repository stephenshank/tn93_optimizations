const { watch } = require('fs');
const { spawn } = require('child_process');


watch('paper/Draft.tex', () => {
  console.log('Building paper at', Date(), '...');


  const snakemake = spawn('snakemake', ['-f', 'paper']);

  snakemake.stderr.on('data', (data) => {
    console.log(`stderr: ${data}`);
  });

  snakemake.on('close', (code) => {
    console.log('... done!');
  });
});

