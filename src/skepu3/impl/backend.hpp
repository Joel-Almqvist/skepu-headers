#pragma once

#include <vector>
#include <iostream>
#include <utility>
#include <cassert>
#include <algorithm>
#include <map>

#include "../backend/debug.h"
#include "skepu3/impl/execution_model.hpp"

#ifndef MAX_DEVICES
#define MAX_DEVICES 4
#endif

#ifndef CPU_THREADS
#define CPU_THREADS 16
#endif

#ifndef GPU_BLOCKS
#define GPU_BLOCKS 256
#endif

#ifndef GPU_THREADS
#define GPU_THREADS 65535
#endif

#ifndef CPU_PARTITION_RATIO
#define CPU_PARTITION_RATIO 0.2f
#endif

namespace skepu
{
	struct Backend
	{
		enum class Type
		{
			Auto, CPU, OpenMP, OpenCL, CUDA, Hybrid
		};
		
		static const std::vector<Type> &allTypes()
		{
			static const std::vector<Backend::Type> types
			{
				Backend::Type::CPU, Backend::Type::OpenMP, Backend::Type::OpenCL, Backend::Type::CUDA, Backend::Type::Hybrid
			};
			
			return types;
		}
		
		static const std::vector<Type> &availableTypes()
		{
			static const std::vector<Backend::Type> types
			{
				Backend::Type::CPU,
#ifdef SKEPU_OPENMP
				Backend::Type::OpenMP,
#endif
#ifdef SKEPU_OPENCL
				Backend::Type::OpenCL,
#endif
#ifdef SKEPU_CUDA
				Backend::Type::CUDA,
#endif
#ifdef SKEPU_HYBRID
                Backend::Type::Hybrid,
#endif
			};
			
			return types;
		}
		
		
		static Type typeFromString(std::string s)
		{
			std::transform(s.begin(), s.end(), s.begin(), ::tolower);
			if (s == "cpu") return Type::CPU;
			else if (s == "openmp") return Type::OpenMP;
			else if (s == "opencl") return Type::OpenCL;
			else if (s == "cuda") return Type::CUDA;
			else if (s == "hybrid") return Type::Hybrid;
			else if (s == "auto") return Type::CUDA;
			else SKEPU_ERROR("Invalid string for backend type conversion");
		}
		
		static bool isTypeAvailable(Type type)
		{
			return type == Backend::Type::Auto ||
				std::find(availableTypes().begin(), availableTypes().end(), type) != availableTypes().end();
		}
		
	};
	
	

	
	
	inline std::ostream &operator<<(std::ostream &o, Backend::Type b)
	{
		switch (b)
		{
		case Backend::Type::CPU:    o << "CPU"; break;
		case Backend::Type::OpenMP: o << "OpenMP"; break;
		case Backend::Type::OpenCL: o << "OpenCL"; break;
		case Backend::Type::CUDA:   o << "CUDA"; break;
		case Backend::Type::Hybrid:   o << "Hybrid"; break;
		case Backend::Type::Auto:   o << "Auto"; break;
		default: o << ("Invalid backend type");
		}
		return o;
	}
	
	
	// "Tagged union" structure for specifying backend and parameters in a skeleton invocation
	struct BackendSpec
	{
		static constexpr Backend::Type defaultType
		{
#if defined(SKEPU_OPENCL)
			Backend::Type::OpenCL
#elif defined(SKEPU_CUDA)
			Backend::Type::CUDA
#elif defined(SKEPU_OPENMP)
			Backend::Type::OpenMP
#else
			Backend::Type::CPU
#endif
		};
		
		static const size_t defaultNumDevices {MAX_DEVICES};
		static const size_t defaultCPUThreads {CPU_THREADS};
		static const size_t defaultGPUThreads {GPU_BLOCKS};
		static const size_t defaultGPUBlocks {GPU_THREADS};
		constexpr static const float defaultCPUPartitionRatio {CPU_PARTITION_RATIO};
		
		BackendSpec(Backend::Type b)
		: m_backend(b) {}
		
		BackendSpec()
		: m_backend(defaultType) {}
		
		size_t devices() const
		{
			return this->m_devices;
		}
		
		void setDevices(size_t numDevices)
		{
			this->m_devices = numDevices;
		}
		
		
		size_t CPUThreads() const
		{
			return this->m_CPUThreads;
		}
		
		void setCPUThreads(size_t threads)
		{
			this->m_CPUThreads = threads;
		}
		
		
		size_t GPUThreads() const
		{
			return this->m_GPUThreads;
		}
		
		void setGPUThreads(size_t threads)
		{
			this->m_GPUThreads = threads;
		}
		
		
		size_t GPUBlocks() const
		{
			return this->m_blocks;
		}
		
		void setGPUBlocks(size_t blocks)
		{
			this->m_blocks = blocks;
		}
		
		
		/**
		 * The partition ratio for hybrid execution, defining how many percent of the workload should be executed by the CPU.
		 */
		float CPUPartitionRatio() const
		{
			return this->m_cpuPartitionRatio;
		}
		
		void setCPUPartitionRatio(float percentage)
		{
			assert(percentage >= 0.0f and percentage <= 1.0f);
			this->m_cpuPartitionRatio = percentage;
		}
		
		Backend::Type backend() const
		{
			return (this->m_backend != Backend::Type::Auto) ? this->m_backend : defaultType;
		}
		
	private:
		Backend::Type m_backend;
		size_t m_devices {defaultNumDevices};
		size_t m_CPUThreads {defaultCPUThreads};
		size_t m_GPUThreads {defaultGPUThreads};
		size_t m_blocks {defaultGPUBlocks};
		float m_cpuPartitionRatio {defaultCPUPartitionRatio};
		
	};
	
	
	/*!
	*  \class ExecPlan
	*
	*  \brief A class that describes an execution plan
	*
	*  This class is used to specifiy execution parameters. For the GPU back ends
	*  you can set both the block size (maxThreads) and the grid size (maxBlocks).
	*  For OpenMP the number of threads is parameterized (numOmpThreads).
	*
	*  It is also possible to specify which backend should be used for a certain
	*  data size. This is done by adding a lowBound and a highBound of data sizes
	*  and a backend that should be used for that range to a list. The skeletons
	*  will use this list when deciding which back end to use.
	*/
	class ExecPlan
	{
	public:
		ExecPlan() : m_calibrated(false)
		{
			this->m_cacheEntry.first = 0;
			cpuModel = nullptr;
			gpuModel = nullptr;
		}
		
		ExecPlan(skepu::ExecutionTimeModel* _cpuModel, skepu::ExecutionTimeModel* _gpuModel) : cpuModel{_cpuModel}, gpuModel{_gpuModel}, m_calibrated(false)
		{
			this->m_cacheEntry.first = 0;
		}
		
		~ExecPlan()
		{
			if(cpuModel)
				delete cpuModel;
			if(gpuModel)
				delete gpuModel;
		}
		
		void add(size_t lowBound, size_t highBound, Backend::Type backend, size_t gs, size_t bs)
		{
			BackendSpec bspec(backend);
			bspec.setGPUThreads(bs);
			bspec.setGPUBlocks(gs);
			this->m_sizePlan.insert(std::make_pair(std::make_pair(lowBound, highBound), bspec));
		}
		
		void add(size_t lowBound, size_t highBound, Backend::Type backend, size_t numOmpThreads)
		{
			BackendSpec bspec(backend);
			bspec.setCPUThreads(numOmpThreads);
			this->m_sizePlan.insert(std::make_pair(std::make_pair(lowBound, highBound), bspec));
		}
		
		void add(size_t lowBound, size_t highBound, BackendSpec bspec)
		{
			this->m_sizePlan.insert(std::make_pair(std::make_pair(lowBound, highBound), bspec));
		}
		
		void add(size_t lowBound, size_t highBound, Backend::Type backend)
		{
			BackendSpec bspec(backend);
			this->m_sizePlan.insert(std::make_pair(std::make_pair(lowBound, highBound), bspec));
		}
		
		
		void setCPUThreads(size_t size, size_t maxthreads)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			this->find(size).setCPUThreads(maxthreads);
		}
		
		size_t CPUThreads(size_t size)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			return this->find(size).CPUThreads();
		}
		
		
		void setGPUThreads(size_t size, size_t maxthreads)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			this->find(size).setGPUThreads(maxthreads);
		}
		
		size_t GPUThreads(size_t size)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			return this->find(size).GPUThreads();
		}
		
		
		void setGPUBlocks(size_t size, size_t maxblocks)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			this->find(size).setGPUBlocks(maxblocks);
		}
		
		size_t GPUblocks(size_t size)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			return this->find(size).GPUBlocks();
		}
		
		
		void setDevices(size_t size, size_t maxdevices)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			this->find(size).setDevices(maxdevices);
		}
		
		size_t devices(size_t size)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			return this->find(size).devices();
		}
		
		
		bool isTrainedFor(size_t size)
		{
			if (this->m_sizePlan.empty())
				return false;
			
			for (auto plan : this->m_sizePlan)
			{
				if (size >= plan.first.first && size <= plan.first.second)
					return true;
			}
			return false;
		}
		
		bool isCalibrated()
		{
			return this->m_calibrated;
		}
		
		void setCalibrated()
		{
			this->m_calibrated = true;
		}
		
		BackendSpec &find(size_t size)
		{
			if (this->m_sizePlan.empty())
				SKEPU_ERROR("Empty execution plan!");
			
			if (this->m_cacheEntry.first == size)
				return m_cacheEntry.second;
			
			for (auto plan : this->m_sizePlan)
			{
				if (size >= plan.first.first && size <= plan.first.second)
				{
					if(cpuModel and gpuModel)
					{
						float cpuRatio = ExecutionTimeModel::predictCPUPartitionRatio(*cpuModel, *gpuModel, size);
						plan.second.setCPUPartitionRatio(cpuRatio);
					}
					
					this->m_cacheEntry = std::make_pair(size, plan.second);
					return this->m_cacheEntry.second;
				}
			}
			
			if(cpuModel and gpuModel)
			{
				float cpuRatio = ExecutionTimeModel::predictCPUPartitionRatio(*cpuModel, *gpuModel, size);
				this->m_sizePlan.rbegin()->second.setCPUPartitionRatio(cpuRatio);
			}
			this->m_cacheEntry = std::make_pair(size, this->m_sizePlan.rbegin()->second);
			return this->m_cacheEntry.second;
		}
		
		void clear()
		{
			this->m_sizePlan.clear();
		}
		
	private:
		std::pair<size_t, BackendSpec> m_cacheEntry;
		std::map<std::pair<size_t, size_t>, BackendSpec> m_sizePlan;
		
		skepu::ExecutionTimeModel* cpuModel;
		skepu::ExecutionTimeModel* gpuModel;
		
		/*! boolean field to specify if this exec plan is properly initialized or not */
		bool m_calibrated;
	};
	
}