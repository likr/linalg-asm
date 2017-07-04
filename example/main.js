/* global fetch, WebAssembly */

const importObject = {
  env: {
    memory: new WebAssembly.Memory({
      initial: 256,
      maximum: 256
    }),
    table: new WebAssembly.Table({
      initial: 0,
      maximum: 0,
      element: 'anyfunc'
    }),
    memoryBase: 0,
    tableBase: 0
  }
}

fetch('./linalg.wasm')
  .then((response) => response.arrayBuffer())
  .then((bytes) => WebAssembly.instantiate(bytes, importObject))
  .then((result) => result.instance)
  .then((linalg) => {
    const n = 100
    const heap = importObject.env.memory.buffer
    const a = new Float64Array(heap, 0, n * n)
    const b = new Float64Array(heap, 8 * n * n, n * n)
    const c = new Float64Array(heap, 16 * n * n, n * n)

    for (let i = 0; i < n; ++i) {
      for (let j = 0; j < n; ++j) {
        a[i * n + j] = Math.random()
        b[i * n + j] = Math.random()
        c[i * n + j] = 0.0
      }
    }

    linalg.exports.dgemm(0, 0, n, n, n, 1.0, a.byteOffset, n, b.byteOffset, n, 0.0, c.byteOffset, n)
    console.log(c)
  })
